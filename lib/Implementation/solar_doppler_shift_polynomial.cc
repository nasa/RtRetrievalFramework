#include "solar_doppler_shift_polynomial.h"
#include "old_constant.h"
#include "wgs84_constant.h"
#include "fp_exception.h"
#include "level_1b.h"
#include "level_1b_fts.h"
#include <vector>

using namespace FullPhysics;
using namespace boost::posix_time;
using namespace boost::gregorian;
using namespace units;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<SolarDopplerShift> create_from_l1b(const Level1b& l1b, int spec_index, bool do_doppler_shift)
{
  return boost::shared_ptr<SolarDopplerShift>(
         new SolarDopplerShiftPolynomial(l1b.time(spec_index),
                                         l1b.latitude(spec_index),
                                         l1b.solar_zenith(spec_index),
                                         l1b.solar_azimuth(spec_index),
                                         l1b.altitude(spec_index),
                                         DefaultConstant(),
                                         do_doppler_shift));
}

boost::shared_ptr<SolarDopplerShift> create_from_runlog(const Level1bFts& l1b, int spec_index, bool do_doppler_shift)
{
  double ds = l1b.run_log(spec_index).observer_sun_doppler_shift * 1e-6;
  return boost::shared_ptr<SolarDopplerShift>
    (new SolarDopplerShiftPolynomial(ds, l1b.time(spec_index), 
                                     DefaultConstant(), do_doppler_shift));
}

REGISTER_LUA_DERIVED_CLASS(SolarDopplerShiftPolynomial, SolarDopplerShift)
.scope
[
 luabind::def("create_from_l1b", &create_from_l1b),
 luabind::def("create_from_runlog", &create_from_runlog)
]
REGISTER_LUA_END()
#endif


//-----------------------------------------------------------------------
/// This is the only place that uses the library CppAD. We have removed
/// the library, but copy this class here to use. This can easily be
/// replaced at some point, but for now I don't want to mess with
/// rewriting and testing the polynomial.
//-----------------------------------------------------------------------

template <class Type, class Vector>
Type Poly(size_t k, const Vector &a, const Type &z)
{        size_t i;
        size_t d = a.size() - 1;

        Type tmp;

        // case where derivative order greater than degree of polynomial
        if( k > d )
        {        tmp = 0;
                return tmp;
        }
        // case where we are evaluating a derivative
        if( k > 0 )
        {        // initialize factor as (k-1) !
                size_t factor = 1;
                for(i = 2; i < k; i++)
                        factor *= i;

                // set b to coefficient vector corresponding to derivative
                Vector b(d - k + 1);
                for(i = k; i <= d; i++)
                {        factor   *= i;
                        tmp       = factor;
                        b[i - k]  = a[i] * tmp; 
                        factor   /= (i - k + 1);
                }
                // value of derivative polynomial
                return Poly(0, b, z);
        }
        // case where we are evaluating the original polynomial
        Type sum = a[d];
        i        = d;
        while(i > 0)
        {        sum *= z;
                sum += a[--i];
        }
        return sum;
}

//-----------------------------------------------------------------------
/// Create a SolarDopplerShiftPolynomial.
///
/// \param t Time of sounding.
/// \param lat The geodetic latitude on the Earth's surface
/// \param sol_zen Solar zenith angle
/// \param sol_az  Solar azmimuth angle (direction looking FOV->sun) in 
///               degrees east of north
/// \param elevation Elevation of point
/// \param constant Instance of Constant subtype, used for
///               solar_distance_param
/// \param apply_doppler_shift Indicates if we should apply the
///               doppler shift or not. If this is false, then we just 
///               pass the wave numbers through unchanged in
///               doppler_stretch. 
//-----------------------------------------------------------------------

SolarDopplerShiftPolynomial::SolarDopplerShiftPolynomial
(
 const Time& t, 
 const DoubleWithUnit& lat, const DoubleWithUnit& sol_zen,
 const DoubleWithUnit& sol_az, const DoubleWithUnit& elevation,
 const Constant& constant,
 bool apply_doppler_shift)
  : apply_doppler_shift_(apply_doppler_shift),
    calculated_doppler_shift(true)
{
  range_check(lat.convert(deg).value, -90.0, 90.0);
  range_check(sol_az.convert(deg).value, 0.0, 360.0);
  range_check(sol_zen.convert(deg).value, -90.0, 90.0);
  calc_solar_distance(constant, t);

  // The "1" here refers to the "1st" derivative.

  const double f = 
    ((OldConstant::wgs84_a - OldConstant::wgs84_b) / OldConstant::wgs84_a).value;
  double slat = sin(lat.convert(rad).value);
  double clat = cos(lat.convert(rad).value);
  double r1 = 1.0 / sqrt(1 - (2*f - f * f) * slat * slat);
  DoubleWithUnit r2 = (OldConstant::wgs84_a * r1 + elevation) * clat / 
    DoubleWithUnit(1, rad);
  doppler_rot_earth_sun_ =
    -1 * OldConstant::earth_rot_speed * r2 * sin(sol_zen.convert(rad).value) * 
    cos(sol_az.convert(rad).value - pi / 2);
  doppler_shift_ = (total_velocity() / 
                    OldConstant::speed_of_light).convert(dimensionless).value;
}

//-----------------------------------------------------------------------
/// Create a SolarDopplerShiftPolynomial. This variation explicitly
/// takes the doppler shift to apply. This is useful for FTS, where we
/// want to use the doppler shift supplied in the runlog rather than
/// calculating our own. Note that the solar distance calculation is
/// still done by the fixed polynomial, only the doppler shift is
/// changed. 
///
/// \param t Time of sounding.
/// \param Doppler_shift The doppler shift factor to use.
/// \param constant Instance of Constant subtype, used for
///               solar_distance_param
/// \param apply_doppler_shift Indicates if we should apply the
///               doppler shift or not. If this is false, then we just 
///               pass the wave numbers through unchanged in
///               doppler_stretch. 
//-----------------------------------------------------------------------

SolarDopplerShiftPolynomial::SolarDopplerShiftPolynomial
(double Doppler_shift,
 const Time& t, 
 const Constant& constant,
 bool apply_doppler_shift)
  : doppler_shift_(Doppler_shift),
    apply_doppler_shift_(apply_doppler_shift),
    calculated_doppler_shift(false)
{
  calc_solar_distance(constant, t);
}

//-----------------------------------------------------------------------
/// Calculate the solar distance
//-----------------------------------------------------------------------

void SolarDopplerShiftPolynomial::calc_solar_distance
(const Constant& constant, const Time& t)
{
  boost::array<double, 7> sdp = constant.solar_distance_param();
  std::vector<double> solar_distance_param_v(sdp.begin(), sdp.end());
  // The polynomial is relative to day of the year, with noon
  // January 1 set to 1. 
  ptime pt(t);
  double doy = pt.date().day_of_year() +
    ((double) pt.time_of_day().total_nanoseconds()) / 
    (hours(24).total_nanoseconds()) + 0.5;
  
  // The "0" here refers to the "0th" derivative, e.g., we just return
  // the results of the polynomial.
  solar_distance_.value = Poly(0, solar_distance_param_v, doy);
  solar_distance_.units = OldConstant::AU;

  solar_velocity_ = DoubleWithUnit(Poly(1, solar_distance_param_v, doy),
                                   OldConstant::AU / day);
}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void SolarDopplerShiftPolynomial::print(std::ostream& Os) const
{
  Os << "SolarDopplerShiftPolynomial\n"
     << "  Solar distance: " << solar_distance_ << "\n"
     << "  Doppler shift:  " << doppler_shift_ << "\n"
     << "  Doppler shift " << (calculated_doppler_shift ? "was calculated" :
                               "was read in") << "\n"
     << "  Apply Doppler:  " << (apply_doppler_shift_ ? "true" : "false");
}

// See base class for description
SpectralDomain SolarDopplerShiftPolynomial::doppler_stretch(
const SpectralDomain& Spec_domain) const
{
  if(apply_doppler_shift_) {
    // The correction is either a multiplication for wavenumbers,
    // or a division for wavelength.
    if(Spec_domain.type_preference() == SpectralDomain::PREFER_WAVENUMBER)
      return SpectralDomain(Array<double, 1>(Spec_domain.data() * 
                                             (1 + doppler_shift_)),
                            Spec_domain.units());
    else
      return SpectralDomain(Array<double, 1>(Spec_domain.data() / 
                                             (1 + doppler_shift_)),
                            Spec_domain.units());
  } else
    return Spec_domain;
}

