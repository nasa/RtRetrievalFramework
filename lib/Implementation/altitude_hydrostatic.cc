#include "altitude_hydrostatic.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AltitudeHydrostatic, Altitude)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
     const boost::shared_ptr<Temperature>&,
     DoubleWithUnit, DoubleWithUnit>())
REGISTER_LUA_END()

#endif


inline double sqr(double x) {return x * x;}
inline AutoDerivative<double> sqr(const AutoDerivative<double>& x) 
{return x * x;}

//-----------------------------------------------------------------------
/// This is a direct transliteration of the old Fortran code to
/// C++. We move this to C++ so we can add the derivative tracking
/// using AutoDerivative.
///
/// From the old Fortran comments:
///  Computes the effective Earth gravity at a given latitude and altitude.
///  This is the sum of the gravitational and centripital accelerations.
///  These are based on equation I.2.4-(17) in US Standard Atmosphere 1962
///  The Earth is assumed to be an oblate ellipsoid, with a ratio of the
///  major to minor axes = sqrt(1+con) where con=.006738
///  This eccentricity makes the Earth's gravititational field smaller at the
///  poles and larger at the equator than if the Earth were a sphere of the
///  same mass. It also makes the local mid-latitude gravity field not point
///  toward the center of mass.
///
///  Input Parameters:
///      gdlat       r*4  GeoDetric Latitude (degrees)
///      altit       r*4  Geometric Altitude (km)
///
///  Output Parameter:
///      gravity     r*4  Effective Gravitational Acceleration (m/s2)
///
///  Interestingly, since the centripital effect of the Earth's rotation
///  (-ve at equator, 0 at poles) has almost the opposite shape to the
///  second order gravitational field (+ve at equator, -ve at poles), their
///  sum is almost constant so that the surface gravity can be approximated
///  (to .07%) by the simple expression g = 0.99746*GM/radius**2, the latitude
///  variation coming entirely from the variation of surface r with latitude.
//-----------------------------------------------------------------------

AutoDerivative<double> AltitudeHydrostatic::gravity_calc(double gdlat, 
                                                         AutoDerivative<double> altit) const
{
  double gclat;                        // geocentric latitude (radians)
  AutoDerivative<double> radius; 
                                // radial distance (metres)
  AutoDerivative<double> ff,hh;
                                // scratch variables

  const double d2r=3.14159265/180;
  const double gm=3.9862216e+14;
  const double omega=7.292116E-05;
  const double con=.006738;
  const double shc=1.6235e-03;
  const double eqrad=6378178;
  const double ge=gm/sqr(eqrad); // = gravity at Re

//  Convert from geodetic latitude (GDLAT) to geocentric latitude (GCLAT).
  gclat=atan(tan(d2r*gdlat)/(1+con));  // radians
  radius=1000*altit+eqrad/sqrt(1.+con*sqr(sin(gclat)));
  ff=sqr(radius/eqrad);
  hh=radius*sqr(omega);
  return (ge*(1-shc*(3*sqr(sin(gclat))-1)/ff)/ff-hh*sqr(cos(gclat)))
    *(1+0.5*sqr(sin(gclat)*cos(gclat)*(hh/ge+2*shc/sqr(ff))));
}

//-----------------------------------------------------------------------
/// Solve the hydrostatic equation to find the altitude and
/// gravitational acceleration at each pressure level.
///
/// Gravity is calculated on a finer grid (sublayers) in order to 
/// reduce create a smoother gravity profile
///
/// Greater accuracy would be obtained by using Tvbar instead of Tbar
/// in the calculation of dz.
/// This would require the passing of a specific humidty profile.
/// Tvbar = Tbar * (1 + eps2 * shbar)
/// Where shbar would be calculated in a similar manner to pbar, tbar
//-----------------------------------------------------------------------

void AltitudeHydrostatic::altitude_calc
(const ArrayAdWithUnit<double, 1>& press_grid, 
 const ArrayAdWithUnit<double, 1>& temp_grid) const
{
  // Gas constant for dry air, in J g^{-1} K^{-1}
  const double Rd = OldConstant::gas_constant / OldConstant::molar_weight_dry_air;

  //const double epsilon = OldConstant::molar_weight_water / OldConstant::molar_weight_dry_air;
  // constant used in calculating virtual temp, ~0.61.
  // Would be used if had access to specific humidity in this routine
  //const double eps2 = (1.0-epsilon)/epsilon;

  Array<AutoDerivative<double>, 1> pg(press_grid.convert(units::Pa).value.to_array());
  ArrayAd<double, 1> tg(temp_grid.convert(units::K).value);

  ArrayAd<double, 1> altv(press_grid.value.rows(), 
                          std::max(press_grid.number_variable(),
                                   temp_grid.number_variable()));

  ArrayAd<double, 1> psublayer((press_grid.rows()-1)*num_sublayer, press_grid.number_variable());
  ArrayAd<double, 1> gravv((altv.rows()-1)*num_sublayer, altv.number_variable());

  // Calculate altitude and gravity at surface using known altitude
  altv(altv.rows() - 1) = (surface_height.value < 0 ? 0 : 
                           surface_height.convert(units::km).value);

  // Calculate altitude and gravity constant
  for(int lev_idx = pg.rows() - 1; lev_idx > 0; lev_idx--) {
    AutoDerivative<double> dp = (pg(lev_idx-1) - pg(lev_idx)) / num_sublayer;
    AutoDerivative<double> dt = (tg(lev_idx-1) - tg(lev_idx)) / num_sublayer;
    AutoDerivative<double> zlo = altv(lev_idx);
    for(int sub_count = 0; sub_count < num_sublayer; sub_count++) {
      int sub_idx = (lev_idx-1) * num_sublayer + num_sublayer - sub_count - 1;
      AutoDerivative<double> plo = pg(lev_idx) + sub_count * dp; 
      psublayer(sub_idx) = plo + 0.5 * dp;
      AutoDerivative<double> logratio = log( plo / ( plo + dp ) );

      AutoDerivative<double> tlo = tg(lev_idx) + sub_count * dt;
      AutoDerivative<double> tbar = tlo + 0.5 * dt;

      gravv(sub_idx) = gravity_calc(latitude.convert(units::deg).value, zlo);
      AutoDerivative<double> dz = (Rd * tbar / gravv(sub_idx)) * logratio;
      gravv(sub_idx) = gravity_calc(latitude.convert(units::deg).value, zlo+0.5*dz);
      zlo = zlo + dz; 
    }
    altv(lev_idx-1) = zlo; 
  }

  std::vector<AutoDerivative<double> > plevlist;
  std::vector<AutoDerivative<double> > altlist;
  for(int i = 0; i < pg.rows(); ++i) {
    plevlist.push_back(pg(i));
    altlist.push_back(altv(i));
  }
  std::vector<AutoDerivative<double> > playlist;
  std::vector<AutoDerivative<double> > glist;
  for(int i = 0; i < psublayer.rows(); ++i) {
    playlist.push_back(psublayer(i));
    glist.push_back(gravv(i));
  }
  alt.reset(new lin_type(plevlist.begin(), plevlist.end(), altlist.begin()));
  grav.reset(new lin_type(playlist.begin(), playlist.end(), glist.begin()));
}

//-----------------------------------------------------------------------
/// Constructor. Latitude for the surface point should be in degrees,
/// and height in meters.
//-----------------------------------------------------------------------

AltitudeHydrostatic::AltitudeHydrostatic
(const boost::shared_ptr<Pressure>& P, const boost::shared_ptr<Temperature>& T,
 const DoubleWithUnit& Latitude, const DoubleWithUnit& Surface_height, const int Num_sublayer)
: cache_is_stale(true),
  latitude(Latitude),
  surface_height(Surface_height),
  p(P),
  t(T),
  num_sublayer(Num_sublayer)
{
  p->add_observer(*this);
  t->add_observer(*this);
}

//-----------------------------------------------------------------------
/// Calculate altitude an gravity by solving hydrostatic
/// equations. These two parameters need to be calculated at the same
/// time. 
///
/// For performance reasons, we cache the results. We only need to
/// calculate if the cache is stale.
//-----------------------------------------------------------------------

void AltitudeHydrostatic::calc_alt_and_grav() const
{
  if(!cache_is_stale)
    return;

  ArrayAdWithUnit<double, 1> pgrid(p->pressure_grid());
  Array<AutoDerivative<double>, 1> tgrid(pgrid.rows());
  for(int i = 0; i < tgrid.rows(); ++i)
    tgrid(i) = t->temperature(pgrid(i)).convert(units::K).value;
  altitude_calc(pgrid, ArrayAdWithUnit<double, 1>(ArrayAd<double, 1>(tgrid),
                                                  units::K));
  cache_is_stale = false;
}

// See base class for description.

boost::shared_ptr<Altitude> AltitudeHydrostatic::clone(
    const boost::shared_ptr<Pressure>& Press,
    const boost::shared_ptr<Temperature>& Temp) const
{
  return boost::shared_ptr<Altitude>(
     new AltitudeHydrostatic(Press, Temp, latitude, surface_height));
}
