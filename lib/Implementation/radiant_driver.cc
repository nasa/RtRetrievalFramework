#include "radiant_driver.h"
#include "linear_algebra.h"
#include "fp_exception.h"
#include "old_constant.h"
#include "wgs84_constant.h"
#include "stokes_coefficient_constant.h"

extern "C"
{
  // arrays  tau,zlev,refrac_lay_grid,tlev,get_atmos_jacobian,l_tau
  // l_user_radiance_direct
  // l_radiance_direct,user_layer,user_tau,l_user_tau,user_radiance_direct

  // in: numlay,fsun,mu0,tau,dim_reset,use_pseudo_spherical new_scene_geo
  //  planet_radius,zlev.use_refraction,refrac_ind_par,refrac_lay_grid
  //  plev,tlev,linearize_atmos_par,get_atmos_jacobian,numpar,l_tau,l_dim_reset
  //  get_user_rad, n_user_rad,user_layer,l_user_tau
  //  user_tau

  // out: radiance_direct,l_radiance_direct,user_radiance_direct,
  // l_user_radiance_direct
  void rad_direct(
		  int *numlay,
		  double *fsun,
		  double *mu0,
		  double *tau,
		  double *radiance_direct, 
		  bool *dim_reset,
		  bool *use_pseudo_spherical,
		  bool *new_scene_geo,
		  double *planet_radius,
		  double *zlev,
		  bool*use_refraction,
		  double *refrac_ind_par,
		  int *refrac_lay_grid,
		  double *plev,
		  double *tlev, 
		  bool *linearize_atmos_par,
		  bool *get_atmos_jacobian,
		  int *numpar,   
		  double *l_tau,
		  bool *l_dim_reset,
		  double *l_radiance_direct,
		  bool *get_user_rad,
		  int *n_user_rad,
		  int *user_layer,
		  double *user_tau,
		  double *l_user_tau,
		  double *user_radiance_direct,
		  double *l_user_radiance_direct);

  void rad_init(double **tau_in_c,double **l_tau_ind, int *maxlay,int *maxatm);
  void rad_cleanup(double *tau_in_c,double *l_tau_ind,int *numpar, int *numlay );

}


//-----------------------------------------------------------------------
/// Constructor.
///
//\param Atm. An atmospsher object.
// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
///////////////////////////////////////

RadiantDriver::RadiantDriver(const boost::shared_ptr<RtAtmosphere>& Atm,
			     const blitz::Array<double, 1>& Sza)
  : sza(Sza.copy())

{
  atm = Atm;
  Array<double, 2> stokes_coef_v(sza.rows(), number_stokes());
  stokes_coef_v = 0;
  stokes_coef_v(Range::all(), 0) = 1.0;
  stokes_coef.reset(new StokesCoefficientConstant(stokes_coef_v));

  for(int i = 0; i < number_spectrometer(); ++i) {
    range_check(sza(i), 0.0, 90.0);
  }

  // Allocates space for radiant arrays
  int maxlay,  maxatm;
  double *tau_ind, *l_tau_ind; 

  rad_init(&tau_ind,&l_tau_ind,&maxlay,&maxatm);
  tau_in.reference(Array<double, 1>(tau_ind, shape(maxlay), neverDeleteData,
				    ColumnMajorArray<1>()));
  
  l_tau_in.reference(Array<double, 2>(l_tau_ind, shape(maxatm, maxlay), 
				      neverDeleteData,
				      ColumnMajorArray<2>()));
  atm->add_observer(*this);
}

//-----------------------------------------------------------------------
/// Destructor
//-----------------------------------------------------------------------

RadiantDriver::~RadiantDriver()
{
  int numpar = l_tau_in.extent(firstDim);
  int numlay = l_tau_in.extent(secondDim);
  rad_cleanup(tau_in.dataFirst(), l_tau_in.dataFirst(), &numpar, &numlay);
}

Array<double, 1> RadiantDriver::stokes_single_wn
(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const {
  return stokes_and_maybe_jacobian(Wn, Spec_index, Iv).value();
}

ArrayAd<double, 1> RadiantDriver::stokes_and_jacobian_single_wn
(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const
{
  return stokes_and_maybe_jacobian(Wn, Spec_index, Iv);
}

ArrayAd<double, 1> RadiantDriver::stokes_and_maybe_jacobian
(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Range rlay(0, atm->number_layer() - 1);
  Range ra(Range::all());
  ArrayAd<double, 1> od, ssa;
  Array<double, 3> jac_iv(0,0,0);

  if(Iv.rows() != 0) {
    od.reference(atm->optical_depth_wrt_iv(Wn, Spec_index, Iv));
    if(!Iv.is_constant())
      jac_iv.reference(Iv.jacobian());
  } else {
    od.reference(atm->optical_depth_wrt_iv(Wn, Spec_index));
    ArrayAd<double, 2> t(atm->intermediate_variable(Wn, Spec_index));
    if(!t.is_constant())
      jac_iv.reference(t.jacobian());
  }
  
  Array<double, 1> jac(jac_iv.depth());
  ArrayAd<double,1> res(number_stokes(), jac.rows());

  tau_in(rlay) = od.value();

  int natm_jac = od.number_variable();
  Range rjac(0, natm_jac - 1);
  l_tau_in(rlay, rjac) = where(od.value()(i1) != 0, 
			       od.jacobian() / od.value()(i1), 0.0);

  int numlay = atm->number_layer(); 
  double fsun = 1.0;
  double mu0 = cos((sza(0)/180.0*OldConstant::pi)); 

  double radiance_direct = 0; // output

  bool dim_reset = true; // this will cause reallocation of memory
  bool use_pseudo_spherical = true; 
  bool new_scene_geo = true; 
  double planet_radius = OldConstant::wgs84_a.convert(units::km).value;

  ArrayAd<double, 1> temp_zlev = atm->altitude(Spec_index).convert(units::km).value;
  blitz::Array<double, 1> zlev = temp_zlev.value();
  bool use_refraction = false; // always false for radiant 

  // matches value from chapman boa currently, if refraction
  // ever enabled
  double refrac_ind_par = 0.000288;

  // hardcoded in old code
  blitz::Array<int, 1> refrac_lay_grid(numlay);
  refrac_lay_grid = 10; //hardcoded in old code 

  // the following 2 are not used because use_refraction is false for radiant
  blitz::Array<double, 1> plev(numlay+1); // not used w/o refraction
  blitz::Array<double, 1> tlev(numlay+1); // not used w/o refraction
  int numpar = natm_jac; //number of rows in jacobian??

  bool linearize_atmos_par = numpar > 0; // true if creating jacobians
  blitz::Array<bool, 2> get_atmos_jacobian(numpar,numlay);
  get_atmos_jacobian = linearize_atmos_par; // true if we need jacobians, i.e. doing retrievals
  bool l_dim_reset = true; //causes reallocation
  blitz::Array<double, 2> l_radiance_direct(numpar,numlay, blitz::ColumnMajorArray<2>());
  
  bool get_user_rad = false; //for radiant
  int n_user_rad = 10; //??

  blitz::Array<int, 1> user_layer(n_user_rad);
  blitz::Array<double, 1> user_tau(n_user_rad);
  blitz::Array<double, 2>  l_user_tau(numpar,n_user_rad);
  blitz::Array<double, 1>  user_radiance_direct(n_user_rad);
  blitz::Array<double, 3>  l_user_radiance_direct(numpar,numlay,n_user_rad);
  rad_direct(&numlay, // n_active_levels -1, from pressure file, might be modified by conner
	     // sounding info pressure file
	     &fsun, //always 1 for radiant
	     &mu0, // layer 1, new_scene_geo
	     tau_in.dataFirst(), //optical depths taug + taur + taua
	     &radiance_direct, //ouput BEER'S LAW COMPUTATION FOR DIRECT RADIANCE
	     &dim_reset, //I'm setting to always true to simplify dome logic
	     &use_pseudo_spherical, //true for radiant
	     &new_scene_geo, //I'm setting to always true to simplify dome logic
	     &planet_radius, //constant from constants module
	     zlev.dataFirst(), //altitudes
	     &use_refraction, //always false for radiant
	     &refrac_ind_par, //hardcoded in old code
	     refrac_lay_grid.dataFirst(), //hardcoded in old code
	     plev.dataFirst(), //not used because no refraction
	     tlev.dataFirst(), // not used because no refraction
	     &linearize_atmos_par,
	     get_atmos_jacobian.data(),
	     &numpar,
	     l_tau_in.dataFirst(),
	     &l_dim_reset,
	     l_radiance_direct.dataFirst(),
	     &get_user_rad,
	     &n_user_rad,
	     user_layer.dataFirst(),
	     user_tau.dataFirst(),
	     l_user_tau.dataFirst(),
	     user_radiance_direct.dataFirst(),
	     l_user_radiance_direct.dataFirst());
    
  radiance_direct *= SOLID_ANGLE;
  l_radiance_direct = l_radiance_direct * SOLID_ANGLE;
  for(int i = 0; i < jac.rows(); ++i) {
    double val = 0;
    for(int m = 0; m < jac_iv.rows(); ++m)
      for(int n = 0; n < jac_iv.cols(); ++n) {
	//note the transposition
	val += l_radiance_direct(n,m) * jac_iv(m, n, i); //this is atmoswf
      }
    jac(i) = val;
  }
  
  res = 0;
  res(0) =  AutoDerivative<double>(radiance_direct,jac);

  return res;
}


//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void RadiantDriver::print(std::ostream& Os) const
{
  Os << "RadiantDriver";
}

