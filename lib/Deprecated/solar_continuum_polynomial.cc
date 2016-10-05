#include "solar_continuum_polynomial.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/foreach.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(SolarContinuumPolynomial, SolarContinuumSpectrum)
.def(luabind::constructor<ArrayWithUnit<double, 1>, bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor. There are two conventions used for the reported solar
/// spectrum. For OCO, we give the results as ph/s/m^2/micron, while
/// for GOSAT we convert using the number of photons at each
/// wavenumber to give W/m^2/cm^-1. You can choose the convention used
/// by specifying if we Convert_from_photon or not.
//-----------------------------------------------------------------------

SolarContinuumPolynomial::SolarContinuumPolynomial(
const ArrayWithUnit<double, 1>& Param,
bool Convert_from_photon)
: param(Param), convert_from_photon(Convert_from_photon)
{
  using namespace units;
  range_min_check((int) Param.value.size(), 1);
  Param.units.is_commensurate(ph / (s * m * m * micron));
}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void SolarContinuumPolynomial::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "    ");
  Os << "SolarContinuumPolynomial:\n"
     << "  Parameters:\n";
  opad << param.value;
  opad.strict_sync();
  Os << "\n  Parameter units: \n";
  opad << param.units;
  opad.strict_sync();
}

// See base class for description.
Spectrum SolarContinuumPolynomial::solar_continuum_spectrum(
  const SpectralDomain& spec_domain) const
{
  // The polynomial is in terms of the wavelength in microns. 
  Array<double, 1> wl = spec_domain.wavelength();

  // Calculate polynomial in wl
  ArrayWithUnit<double, 1> res;
  res.value.resize(wl.shape());
  res.units = param.units;
  res.value = param.value(param.value.size() - 1);
  for(int i = ((int) param.value.size()) - 2; i >= 0; i--) {
    res.value *= wl;
    res.value += param.value(i);
  }

  // Convert from ph/s/m^2/micron to W/m^2/cm^-1 if requested
  if(convert_from_photon)
    res *= spec_domain.photon_to_radiance_factor();
  return Spectrum(spec_domain, SpectralRange(res.value, res.units));
}
