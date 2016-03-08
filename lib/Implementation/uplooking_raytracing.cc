#include "uplooking_raytracing.h"
#include "old_constant.h"
#include "wgs84_constant.h"
#include "array_with_unit.h"
#include "logger.h"

using namespace FullPhysics;
using namespace blitz;

Array<double, 1> calculate_apparent_surface_sza_value
(const boost::shared_ptr<Level1b>& L1b,
 const boost::shared_ptr<RtAtmosphere>& Atm)
{
  const boost::shared_ptr<Level1bFts> l1b_fts(boost::dynamic_pointer_cast<Level1bFts>(L1b));
  const boost::shared_ptr<AtmosphereOco> atm_oco(boost::dynamic_pointer_cast<AtmosphereOco>(Atm));
  return UplookingRaytracing::calculate_apparent_surface_sza(l1b_fts, atm_oco).value;
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(UplookingRaytracing)
.scope
[
 luabind::def("calculate_apparent_surface_sza", &calculate_apparent_surface_sza_value)
]
REGISTER_LUA_END()
#endif

ArrayWithUnit<double, 1> UplookingRaytracing::calculate_apparent_surface_sza
(const boost::shared_ptr<Level1bFts>& l1b,
 const boost::shared_ptr<AtmosphereOco>& atm)
{
  ArrayWithUnit<double, 1> sza;
  sza.value.resize(l1b->number_spectrometer());
  sza.value = 0.0;
  sza.units = units::deg;

  Array<AutoDerivative<double>, 1> press_grid( atm->pressure_ptr()->pressure_grid().convert(units::Pa).value.to_array() );
  Array<AutoDerivative<double>, 1> temp_grid(press_grid.rows());
  for(int i = 0; i < temp_grid.rows(); ++i)
    temp_grid(i) = 
      atm->temperature_ptr()->temperature(AutoDerivativeWithUnit<double>(press_grid(i), units::Pa)).convert(units::K).value;
  int const nlayers = temp_grid.extent(firstDim)-1;

  for(int Ispec = 0; Ispec < l1b->number_spectrometer(); ++Ispec) {
    if (l1b->radiance(Ispec).data().rows() > 0) {
      Array<AutoDerivative<double>, 1> height_grid( atm->altitude(Ispec).convert(units::km).value.to_array());

      FtsRunLogRecord fr = l1b->run_log(Ispec);

      double asza_in = fr.solar_zenith + fr.zenith_offset;
      double fovr = 90.0/3.14159265 * fr.internal_fov;
      double wavtkr = fr.sun_tracker_frequency;
      double wavmic = 0.5*fr.spacing_raw_spectrum*(fr.index_first + fr.index_last);
      double Psurf = fr.outside_pressure / 1013.25; // in atm

      double * z_loc = new double[nlayers+1];
      double * t_loc = new double[nlayers+1];
      double * p_loc = new double[nlayers+1];
      double * sp = new double[nlayers+1];

      for (int i=0; i<=nlayers; i++) {
	t_loc[nlayers-i] = temp_grid(i).value();
	z_loc[nlayers-i] = height_grid(i).value();
	p_loc[nlayers-i] = press_grid(i).value() / 101312.5;  // Pa to atm
      }

      // find height (zobs): i_lay is first level above surface
      int i_lay=0;
      while(i_lay<=nlayers) {
	if (p_loc[i_lay] <= Psurf) break;
	i_lay++;
      }

      double zobs = 0.0;
      if (i_lay > nlayers)  {
	Logger::info() << "UplookingRaytracing: Warning: all levels are below surface!\n";
	i_lay=nlayers;
	zobs=z_loc[i_lay];
      }
      else {
	double fr = log(p_loc[i_lay-1]/Psurf)/log(p_loc[i_lay-1]/p_loc[i_lay]);
	double fta=1.0;
	double ftawas;

	int it=0;
	for (; it<10; it++) {
	  ftawas=fta;
	  double tobs=t_loc[i_lay-1]+fta*fr*(t_loc[i_lay]-t_loc[i_lay-1]);
	  fta=(t_loc[i_lay-1]+tobs)/(t_loc[i_lay-1]+t_loc[i_lay]);
	  if(abs(fta-ftawas) < 2.e-7)  break;
	}
	if (it==10)
	  Logger::info() << "UplookingRaytracing: WARNING: unconverged\n";
	zobs=z_loc[i_lay-1]+fta*fr*(z_loc[i_lay]-z_loc[i_lay-1]);
      }

      double roc = OldConstant::wgs84_a.convert(units::km).value;

      int ifail;
      double zmin, bend;
      tlpath(nlayers+1, z_loc, t_loc, p_loc, asza_in, fovr, roc, zobs, wavtkr, wavmic, &zmin, &bend, sp, &ifail);

      sza.value(Ispec) = l1b->solar_zenith(Ispec).convert(units::deg).value - bend;

      delete [] t_loc;
      delete [] p_loc;
      delete [] z_loc;
    }
  }

  return sza;
}

inline double calc_opcon(double wav) {
  //  Precompute optical constant used to calculate refractive index gradient
  double xt=wav*wav*1.0e-8;
  double a=64.328e0;
  double b=29498.1e0/(146.0e0-xt);
  double c=255.4e0/(41.0e0-xt);
  double aa=1.0e-6*(a+b+c);
   
  //return 0.0;  //  uncomment to disable refraction
  return aa*(aa+2.0e0)/(1.225014e0*(aa*(aa+2.0e0)+3.0e0));
}

inline double ri(double temp, double pres, double opcon) {

  double d=353.0*pres*opcon/temp;

  // cut off refractive index
  return (d<=0.97e0 ? sqrt((1.0+2.0*d)/(1.0-d)) : 9.9);
}

inline double grad(double temp, double pres, double rpsh, double opcon) {
  //  Calculates  1/n.dn/dz  where n=refractive index.
  //  rpsh = reciprocal scale height
  double d=353.0*pres*opcon/temp;

  return 1.5*d*rpsh/((1.0+2.0*d)*(1.0-d));
}

void find(int nlev, double * z, double zhp, int * lev1, int * lev2) {
  //  Finds the model levels which bracket the altitude zhp.
  int l2;

  for (l2=0; l2<nlev; ++l2) {
    if(z[l2]>=zhp) break;
  }

  if (l2==nlev) {
    cout << "WARNING: ray above TOA" << endl;
    --l2;
  }

  if (l2==0) {
    cout << "WARNING: ray below BOA" << endl;
    ++l2;
  }

  *lev1 = l2-1;
  *lev2 = l2;
}

void getsp(int nlev, double * z, double * p, double * t, double ths, double roc,
	   double zobs, double * zmin, double opcon, double * spg,
	   double * bend, int * ifail)
{
  //  Calculates the geometrical slant paths sp(nlev) for a ray of apparent
  //  zenith angle ths degrees, hitting an observer at an radius zobs km,
  //  from the centre of curvature.
  //
  //  The atmospheric temperature and pressure t(nlev) & p(nlev) is tabulated
  //  at radii z(nlev), which need not be equally spaced.
  //
  //  The subroutine also returns bend, the angle (degrees) through which the
  //  ray was deviated by refraction. this permits getsp to be put into an
  //  iteration loop tp find the angle ths which gives rise to a particular
  //  astronomical solar zenith angle.

  int const maxstep=64000;
  int lev1,lev2,k;
  double ztoa,zhp,del_zlev,pres,temp,x1,x2,a1,a2,con;
  double ds,dx,dz,dt,dphi,phi,th,thp,zh,rpsh,rp;

  const double pi=OldConstant::pi;
  const double d2r=FullPhysics::conversion(units::deg, units::rad);
  const double piby2=0.5*pi;
   
  for (k=0; k<nlev; ++k)
    spg[k]=0.0;

  *zmin=zobs;
  *bend=0.0;
  ztoa=z[nlev-1];

  //  If observer is outside atmosphere, start integration at top of atmosphere
  //  assuming no refraction above Ztoa
  th=ths*d2r;
  if(zobs>ztoa) {
    if(th<piby2) 
      return;
    zh=double(ztoa);
    *zmin=roc*(sin(th)-1.0e0) + sin(th)*zobs;
    *zmin=sin(th)*(roc+zobs)-roc;

    if( *zmin>ztoa ) {
      *ifail=2;  //  ray did not enter atmosphere
      return;
    }

    th=pi-asin(sin(th)*(roc+zobs)/(roc+zh));
    phi=ths*d2r-th;
  } 
  else {
    zh=double(zobs);
    phi=0.0e0;
  }

  //  ds = integration step  < (layer thickness)/8, or 0.5 km horizontal.
  find(nlev,z,zh,&lev1,&lev2);

  dx=0.5;                    // max allowed horizontal step 
  //      dx=sqrt(2.5*ztoa*roc)/maxstep // max allowed horizontal step 
  dz=(z[lev1]-z[lev2])/4;    // max allowed vertical step 
  ds=1.0e0 / (1/dz>abs(cos(th)/dz) ? 1/dz : abs(cos(th)/dz));
  dz=ds*cos(th);
  dx=ds*sin(th);

  zhp=zh+dz/2.;
  dphi=dx/(roc+zhp);
  find(nlev,z,zhp,&lev1,&lev2);
  del_zlev=z[lev1]-z[lev2];
  temp=t[lev1] - (z[lev1]-zhp)*(t[lev1]-t[lev2])/del_zlev;
  rpsh=log(p[lev2]/p[lev1])/del_zlev;
  pres=p[lev1]*exp((z[lev1]-zhp)*rpsh);
  dt=dx*grad(temp,pres,rpsh,opcon)-dphi;
  thp=th+dt/2.;

  //  zhp & thp are the altitude and local zenith angle estimated for the middle
  //  of the next path segment.
  for(k=0; k<maxstep; ++k) {
    dz=ds*cos(thp);
    dx=ds*sin(thp);
    dphi=dx/(roc+zhp);
    del_zlev=z[lev1]-z[lev2];
    temp=t[lev1] - (z[lev1]-zhp)*(t[lev1]-t[lev2])/del_zlev;

    //  Pressure is expressed as a linear combination of pressures at levels lev1 & lev2.
    //  the coefficients a1 & a2 are exponential functions of height, so a1+a2 < 1.
    //  a1 & a2 are also the contributions of the current path segment to the slant
    //  paths of levels lev1 & lev2 respectively.
    x1=(z[lev1]-zhp)/del_zlev;
    x2=x1-1.0;
    rp=log(p[lev2]/p[lev1]);
    a1=-(x2*exp(x1*rp));
    a2= x1*exp(x2*rp);
    pres = p[lev2]*a2 + p[lev1]*a1;
    dt=dx*grad(temp,pres,rp/del_zlev,opcon)-dphi;

    spg[lev1]=spg[lev1]+a1*ds;
    spg[lev2]=spg[lev2]+a2*ds;
    phi=phi+dphi;
    zh=zh+dz;
    th=th+dt;

    //  Terminate integration cleanly if a tangent point is detected.
    if( ths>90.0e0 ) {   // if the initial angle was > 90
      if( (th-piby2)*(th-piby2-dt)<=0.0e0) {  // Gone past TP
	// Back-track to tangent point
	con=(th-piby2)/dt;
	spg[lev1]=spg[lev1]-a1*ds*con;
	spg[lev2]=spg[lev2]-a2*ds*con;
	zh=zh-dz*con;
	th=th-dt*con;
	phi=phi-dphi*con;
	break;
      }
    }

    if(zh>z[lev2]) {
      con=(zh-z[lev2])/dz;
      spg[lev1]=spg[lev1]-a1*ds*con;
      spg[lev2]=spg[lev2]-a2*ds*con;
      zh=z[lev2];
      th=th-dt*con;
      phi=phi-dphi*con;
      if(zh>=ztoa)
	break;
    }
    find(nlev,z,zhp,&lev1,&lev2);
    zhp=zh+dz/2.;
    thp=th+dt/2.;
  } // k=1,maxstep

  if (k==maxstep) {
    *ifail=3;
    return;
  }

  *bend=(th+phi)/d2r-ths;
  *zmin=(zobs<zh ? zobs : zh);
}

void UplookingRaytracing::tlpath(int nlev,double * z, double * t, double * p, double asza_in,
				 double fovr, double roc, double zobs, double wavtkr,
				 double wavmic, double * zmin_out, double * bend_out,
				 double * sp, int * ifail)
{
  double & zmin = *zmin_out;
  double & bend = *bend_out;

  int const mlev = 1000;
  int lev1,lev2,it,iter,k,l;
  double spg[mlev];
  double asza,zmin0,bend0,bendx,rsza,del_zlev,
    pobs,tobs,con,ztan,ptan,ttan,prdiff,
    del_ztan,del_rsza,dum,xsza,weight,rpsh,xx,dbdt,b;

  double opcon_tkr,opcon_mic,rizobs,riztan,rtnt;
  double const pi=3.1415926536;
  double const d2r=pi/180.0;
  int const maxiter=mlev/2;

  *ifail=0;
  if(nlev==1) {
    cout << "TLPATH exiting due to NLEV<=1" << endl;
    sp[0]=0.0;
    bend=0.0;
    zmin=z[0];
    *ifail=6;
    return;
  }

  opcon_tkr=calc_opcon(wavtkr);
  rsza=abs(asza_in);
  if(nlev>mlev) {
    cout << "nlev>mlev" << endl;
    exit(1);
  }

  if(asza_in>=0.0) {   // angle is not refracted, iterate
    //  estimate, using a simple empirical model that involves no ray tracing,
    //  the refracted angle, rsza, that gives rise to the astronomical angle = asza
    //  the purpose of this is to reduce the number of iterations that the ray
    //  has to be traced through the atmosphere. this is achieved not just by a
    //  more accurate starting guess, but also better partial differentials, dbdt.
    //  tests have shown that the rsza estimate is good to better than 0.02 degrees
    //  with any model for both ground-based and atmos geometries.
    asza=rsza;
    find(nlev,z,zobs,&lev1,&lev2);

    del_zlev=z[lev1]-z[lev2];
    rpsh=log(p[lev2]/p[lev1])/del_zlev;
    pobs=p[lev1]*exp((z[lev1]-zobs)*rpsh);
    tobs=t[lev1] - (z[lev1]-zobs)*(t[lev1]-t[lev2])/del_zlev;
    rizobs=ri(tobs,pobs,opcon_tkr);
    con=.035;

    for (it=0; it<maxiter; ++it) {
      xx=1.0/(con+abs(cos(rsza*d2r)));
      bend=194.0*pobs*sin(rsza*d2r)*con*xx/tobs;
      dbdt=1.25*bend*xx;
      if(rsza>90.0) {
	// use snell's law in spherical geometry to evaluate true tangent ht.
	rtnt=rizobs*(roc+zobs)*sin(rsza*d2r);
	ztan=rtnt-roc;
	if(ztan<z[1]-99.999) {
	  *ifail=4;
	  return;
	}
	for (int jt=0; jt<maxiter; ++jt) {
	  find(nlev,z,ztan,&lev1,&lev2);
	  del_zlev=z[lev1]-z[lev2];
	  rpsh=log(p[lev2]/p[lev1])/del_zlev;
	  ptan=p[lev1]*exp((z[lev1]-ztan)*rpsh);
	  ttan = t[lev1] - (z[lev1]-ztan)*(t[lev1]-t[lev2])/del_zlev;
	  riztan=ri(ttan,ptan,opcon_tkr);
	  prdiff=riztan+(roc+ztan)*(1.0-riztan)*rpsh;
	  if(abs(prdiff)<0.66)
	    break;
	  del_ztan = (rtnt-(roc+ztan)*riztan)/prdiff;
	  ztan=ztan+del_ztan;
	  if(abs(del_ztan)<0.001)
	    break;
	}
	// assumes ray bending is directly proportional to tangent density.
	b=388.0*ptan/ttan;
	bend=-bend+b;
	dbdt=-dbdt+1.72*b*(roc+zobs)*rpsh/xx;
      }

      del_rsza=asza-rsza-bend;
      if(rsza>25.0) 
	del_rsza=rsza*del_rsza/(rsza+dbdt);
      rsza=rsza+del_rsza;
      if(abs(del_rsza)<0.0003)
	break;
    }

    if (it==maxiter)
      *ifail=5;

    //  iterate to find the apparent (refracted) angle, rsza, such that rays
    //  of wavenumber wavtkr, when traced back through the atmosphere, emerge
    //  at the desired astronomical angle asza. the estimate of the previous
    //  program section is used as a starting guess for the iteration. speed of
    //  convergence is optimised by using partial differentials dbdt calculated
    //  in the previous program section.

    for (iter=0; iter<2*maxiter; ++iter) {
      getsp(nlev,z,p,t,rsza,roc,zobs,&zmin0,opcon_tkr,spg,&bend,ifail);
      if(rsza>90.0) {
	getsp(nlev,z,p,t,180.0-rsza,roc,zobs,&dum,opcon_tkr,spg,&bendx,ifail);
	bend=2.0*bend+bendx;
      }
      del_rsza=asza-rsza-bend;
      if(rsza>25.0) 
	del_rsza=rsza*del_rsza/(rsza+dbdt);
      rsza=rsza+del_rsza;
      if(abs(del_rsza) < 0.00004)
	break;
    }

    if (iter==2*maxiter)
      *ifail=1;
  }

  for (k=0; k<nlev; ++k) {
    sp[k]=0.0;
  }

  //  having found rsza we now trace rays of wavenumber wavmic from the top, middle
  //  and bottom of the field of view and weight their slant contributions by
  //  0.125, 0.750 & 0.125 respectively. this corrects exactly for quadratic
  //  curvature of the atmospheric airmass as a function of angle over the
  //  circular field of view.
  opcon_mic=calc_opcon(wavmic);
  for (l=-1; l<=1; ++l) {
    weight=0.750-0.625*l*l;
    xsza=rsza+l*fovr;
    zmin=zmin0;     // we want the value when l=0, not l=1
    bend=bend0;

    getsp(nlev,z,p,t,xsza,roc,zobs,&zmin0,opcon_mic,spg,&bend0,ifail);

    if(xsza>90.0) {
      for (k=0; k<nlev; ++k) {
	sp[k]=sp[k]+2*weight*spg[k];
      }
      getsp(nlev,z,p,t,180.0-xsza,roc,zobs,&dum,opcon_mic,spg,&bendx,ifail);
      bend0=2*bend0+bendx;
    }
    for (k=0; k<nlev; ++k) {
      sp[k]=sp[k]+weight*spg[k];
    }
  } // l=-1,1

}
