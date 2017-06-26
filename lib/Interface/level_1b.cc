#include "level_1b.h"
#include "fe_disable_exception.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA

// Convenience function to give us solar zenith, relative azimuth, 
// stokes_coefficents, and sounding azimuth in different forms for Lua

blitz::Array<double, 1> level_1b_sza(const Level1b& Lev1)
{
  blitz::Array<double, 1> res(Lev1.number_spectrometer());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = Lev1.solar_zenith(i).convert(units::deg).value;
  return res;
}

double level_1b_sza_i(const Level1b& Lev1, int Spec_index)
{
  return Lev1.solar_zenith(Spec_index).convert(units::deg).value;
}

blitz::Array<double, 1> level_1b_zen(const Level1b& Lev1)
{
  blitz::Array<double, 1> res(Lev1.number_spectrometer());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = Lev1.sounding_zenith(i).convert(units::deg).value;
  return res;
}

double level_1b_zen_i(const Level1b& Lev1, int Spec_index)
{
  return Lev1.sounding_zenith(Spec_index).convert(units::deg).value;
}

blitz::Array<double, 1> level_1b_azm(const Level1b& Lev1)
{
  blitz::Array<double, 1> res(Lev1.number_spectrometer());
  for(int i = 0; i < res.rows(); ++i) {
    res(i) = Lev1.relative_azimuth(i).convert(units::deg).value;
  }
  return res;
}

double level_1b_azm_i(const Level1b& Lev1, int Spec_index)
{
  return Lev1.relative_azimuth(Spec_index).convert(units::deg).value;
}

blitz::Array<double, 2> level_1b_stokes(const Level1b& Lev1)
{
  blitz::Array<double, 2> res(Lev1.number_spectrometer(), 4);
  for(int i = 0; i < res.rows(); ++i)
    res(i, blitz::Range::all()) = Lev1.stokes_coefficient(i);
  return res;
}

blitz::Array<double, 2> level_1b_s_coeffs(const Level1b& Lev1)
{
  int nparam = Lev1.spectral_coefficient(0).value.rows();
  blitz::Array<double, 2> res(Lev1.number_spectrometer(), nparam);
  for(int i = 0; i < res.rows(); ++i)
    res(i, blitz::Range::all()) = Lev1.spectral_coefficient(i).value;
  return res;
}

ArrayWithUnit<double, 2> level_1b_s_coeffs_with_unit(const Level1b& Lev1)
{
  int nparam = Lev1.spectral_coefficient(0).value.rows();
  Unit units = Lev1.spectral_coefficient(0).units;
  blitz::Array<double, 2> res(Lev1.number_spectrometer(), nparam);
  for(int i = 0; i < res.rows(); ++i)
    res(i, blitz::Range::all()) = Lev1.spectral_coefficient(i).value;
  return ArrayWithUnit<double,2>(res, units);
}

blitz::Array<double, 1> level_1b_radiance(const Level1b& Lev1, int Spec_index)
{
  return Lev1.radiance(Spec_index).data();
}

blitz::Array<double, 1> level_1b_uncertainty(const Level1b& Lev1, int Spec_index)
{
  return Lev1.radiance(Spec_index).uncertainty();
}

ArrayWithUnit<double, 1> level_1b_uncertainty_with_unit(const Level1b& Lev1, int Spec_index)
{
  return ArrayWithUnit<double, 1>(Lev1.radiance(Spec_index).uncertainty(), Lev1.radiance(Spec_index).units());
}

ArrayWithUnit<double, 1> level_1b_radiance_with_unit(const Level1b& Lev1, int Spec_index)
{
  return ArrayWithUnit<double, 1>(Lev1.radiance(Spec_index).data(), Lev1.radiance(Spec_index).units());
}

SpectralRange level_1b_radiance_spectral_range(const Level1b& Lev1, int Spec_index)
{
  return Lev1.radiance(Spec_index);
}

double level_1b_relative_velocity(const Level1b& Lev1, int i)
{
  return Lev1.relative_velocity(i).value;
}

int level_1b_year(const Level1b& Lev1, int i)
{
  boost::posix_time::ptime t = Lev1.time(i);
  // +1900 since value is years since 1900
  return 1900 + to_tm(t).tm_year;
}

int level_1b_month(const Level1b& Lev1, int i)
{
  boost::posix_time::ptime t = Lev1.time(i);
  // +1 because tm_mon goes from 0 to 11. We want 1 to 12.
  return to_tm(t).tm_mon + 1;
}

int level_1b_day(const Level1b& Lev1, int i)
{
  boost::posix_time::ptime t = Lev1.time(i);
  return to_tm(t).tm_mday;
}

int level_1b_hour(const Level1b& Lev1, int i)
{
  boost::posix_time::ptime t = Lev1.time(i);
  return to_tm(t).tm_hour;
}

int level_1b_minute(const Level1b& Lev1, int i)
{
  boost::posix_time::ptime t = Lev1.time(i);
  return to_tm(t).tm_min;
}

int level_1b_second(const Level1b& Lev1, int i)
{
  boost::posix_time::ptime t = Lev1.time(i);
  return to_tm(t).tm_sec;
}

int level_1b_dayofyear(const Level1b& Lev1, int i)
{
  boost::posix_time::ptime t = Lev1.time(i);
  return to_tm(t).tm_yday;
}

DoubleWithUnit l1b_signal_no_indexes(const boost::shared_ptr<Level1b>& l1b, int Spec_index)
{
    return l1b->signal(Spec_index, std::vector<int>());
}

#include "register_lua.h"
REGISTER_LUA_CLASS(Level1b)
.def("number_spectrometer", &Level1b::number_spectrometer)
.def("latitude", &Level1b::latitude)
.def("longitude", &Level1b::longitude)
.def("altitude", &Level1b::altitude)
.def("sza", &level_1b_sza)
.def("sza", &level_1b_sza_i)
.def("sza_with_unit", &Level1b::solar_zenith)
.def("azm", &level_1b_azm)
.def("azm", &level_1b_azm_i)
.def("azm_with_unit", &Level1b::relative_azimuth)
.def("zen", &level_1b_zen)
.def("zen", &level_1b_zen_i)
.def("zen_with_unit", &Level1b::sounding_zenith)
.def("stokes_coef", &level_1b_stokes)
.def("spectral_coefficient", &level_1b_s_coeffs)
.def("spectral_coefficient_with_unit", &Level1b::spectral_coefficient)
.def("spectral_coefficient_with_unit", &level_1b_s_coeffs_with_unit)
.def("year", &level_1b_year)
.def("month", &level_1b_month)
.def("day", &level_1b_day)
.def("hour", &level_1b_hour)
.def("minute", &level_1b_minute)
.def("second", &level_1b_second)
.def("dayofyear", &level_1b_dayofyear)
.def("radiance", &level_1b_radiance)
.def("uncertainty", &level_1b_uncertainty)
.def("uncertainty_with_unit", &level_1b_uncertainty_with_unit)
.def("radiance_with_unit", &level_1b_radiance_with_unit)
.def("radiance_spectral_range", &level_1b_radiance_spectral_range)
.def("signal", &Level1b::signal)
.def("signal", &l1b_signal_no_indexes)
.def("relative_velocity", &level_1b_relative_velocity)
.def("time", &Level1b::time)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<Level1b> >::*pbt1)(
        const std::vector<boost::shared_ptr<Level1b> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<Level1b> >,
                        VectorLevel1b)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<Level1b> >::push_back))
REGISTER_LUA_END()

// Azimuth is modified because the convention used by the OCO L1B file is to
// take both solar and observation angles as viewed from an observer
// standing in the FOV.  In this convention, the convention for glint
// would be a relative azimuth difference of 180 degrees, as the
// spacecraft and sun would be on opposite sides of the sky. However, the
// radiative transfer convention is that the azimuth angles must be the
// same for glint (it is "follow the photons" convention). However, we'd
// like the solar azimuth to not be changed, so as to continue to agree
// with zenith, so this change of the observation azimuth has the effect
// of putting everything in a "reverse follow-the-photons" convention,
// where we look from the satellite to the FOV, then from the FOV to the
// sun.  Note that because of an old historical reason, however, both
// zenith angles remain > 0 and < 90, even in the RT convention.

DoubleWithUnit Level1b::relative_azimuth(int Spec_index) const
{
    Unit orig_unit = sounding_azimuth(Spec_index).units;
    double res = (180 + sounding_azimuth(Spec_index).convert(units::deg).value) - 
        solar_azimuth(Spec_index).convert(units::deg).value;

    if (res > 360)
        res -= 360;
    else if(res < 0)
        res += 360;
      
    return DoubleWithUnit(res, units::deg).convert(orig_unit);
}

DoubleWithUnit Level1b::signal(int Spec_index, const std::vector<int>& Sample_indexes) const
{
    // Basically a copy of what is in ErrorAnalysis::signal
    FeDisableException disable_fp;
    const int nrad = 10;
    SpectralRange rad(radiance(Spec_index));
  
    if(rad.data().rows() ==0) {
        return DoubleWithUnit(0, rad.units());
    }
  
    // Copy either the whole data if no sample index list is used
    // or only selected samples indexes listed in Sample_indexes
    Array<double, 1> used_rad;
    if (Sample_indexes.size() == 0) {
        used_rad.resize(rad.data().rows());
        used_rad = rad.data();
    } else {
        used_rad.resize(Sample_indexes.size());
        for(int samp_idx = 0; samp_idx < (int) Sample_indexes.size(); samp_idx++) {
            used_rad(samp_idx) = rad.data()(Sample_indexes[samp_idx]);
        }
    }
    std::sort(used_rad.data(), used_rad.data() + used_rad.rows()); // Min to max value
    used_rad.reverseSelf(firstDim);     // Now max to min value
    Range avg_range(0, std::min(nrad - 1, used_rad.rows() - 1));
    return DoubleWithUnit(sum(used_rad(avg_range) / avg_range.length()), rad.units());
}

#endif
