#include "solar_model.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(SolarModel, SpectrumEffect)
.def("solar_spectrum", &SolarModel::solar_spectrum)
REGISTER_LUA_END()

#endif

//-----------------------------------------------------------------------
/// Apply the solar model.
///
/// \param Spec Spectrum without solar model applied. This is per
///          solid angle, e.g., sr^-1 (RT code generates "sun-normalized" output)
/// \return Spectrum with solar model applied. 
//-----------------------------------------------------------------------

Spectrum SolarModel::apply_solar_model(const Spectrum& Spec) const
{ 
  if(!Spec.spectral_range().units().is_commensurate(units::inv_sr)) {
    Exception e;
    e << "SolarModel::apply_solar_model expects the input spectrum to be in units of\n"
      << "inverse solid angle (e.g., sr^-1). Instead, it was passed a spectrum with\n"
      << "units of:\n"
      << Spec.spectral_range().units() << "\n";
    throw e;
  }
  ArrayAd<double,1> res(Spec.spectral_range().data_ad().copy());
  Spectrum s(solar_spectrum(Spec.spectral_domain()));
  res.value() *= s.spectral_range().data();
  for(int i = 0; i < res.number_variable(); ++i)
    res.jacobian()(blitz::Range::all(), i) *= s.spectral_range().data();
  return Spectrum(Spec.spectral_domain(), SpectralRange(res, 
	s.spectral_range().units() * Spec.spectral_range().units()));
}

