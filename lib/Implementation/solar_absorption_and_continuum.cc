#include "solar_absorption_and_continuum.h"
#include "old_constant.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(SolarAbsorptionAndContinuum, SpectrumEffect)
.def(luabind::constructor<const boost::shared_ptr<SolarDopplerShift>&,
     const boost::shared_ptr<SolarAbsorptionSpectrum>&,
     const boost::shared_ptr<SolarContinuumSpectrum>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void SolarAbsorptionAndContinuum::print(std::ostream& Os) const
{
  Os << "SolarAbsorptionAndContinuum:\n";
  OstreamPad opad(Os, "  ");
  opad << doppler_shift() << "\n"
       << absorption_spectrum() << "\n"
       << continuum_spectrum();
}

// See base class for  description

Spectrum SolarAbsorptionAndContinuum::solar_spectrum(const SpectralDomain& Spec_domain) const
{
  SpectralDomain sd(doppler_shift().doppler_stretch(Spec_domain));
  double solar_dist = 
    doppler_shift().solar_distance().convert(OldConstant::AU).value;
  Spectrum solar_cont(continuum_spectrum().solar_continuum_spectrum(sd)); 
  Spectrum solar_abs(absorption_spectrum().solar_absorption_spectrum(sd));
  if(!solar_abs.spectral_range().units().is_commensurate(units::dimensionless)) {
    Exception e;
    e << "Absorption spectrum should be dimensionless. Instead, it had the\n"
      << "dimensions of:\n" << solar_abs.spectral_range().units();
    throw e;
  }
  Array<double, 1> res(solar_abs.spectral_range().data() *
		       solar_cont.spectral_range().data() / (solar_dist * solar_dist));
  return Spectrum(Spec_domain, 
		  SpectralRange(res, solar_cont.spectral_range().units()));
}
