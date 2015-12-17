#include "instrument_doppler.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(InstrumentDoppler, SpectrumEffect)
.def(luabind::constructor<const DoubleWithUnit&, const bool>())
.def(luabind::constructor<const DoubleWithUnit&>())
.def(luabind::constructor<const double, const std::string&, const bool>())
.def(luabind::constructor<const double, const std::string&>())
REGISTER_LUA_END()
#endif
 
InstrumentDoppler::InstrumentDoppler(const DoubleWithUnit& Relative_velocity, 
                                     const bool Used_flag)
  : vel_units(Relative_velocity.units)
{
  Array<double, 1> c(1);
  c = Relative_velocity.value;
  Array<bool, 1> f(1);
  f = Used_flag;
  init(c, f);
}

InstrumentDoppler::InstrumentDoppler(const double Relative_velocity_value,
                                     const std::string& Relative_velocity_units,
                                     const bool Used_flag) 
  : vel_units(Relative_velocity_units)
{
  Array<double, 1> c(1);
  c = Relative_velocity_value;
  Array<bool, 1> f(1);
  f = Used_flag;
  init(c, f);
}

void InstrumentDoppler::apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const
{
  // Compute the wavelengths of the object in (2) in the spacecraft frame of reference.
  // Let v_ei be the relative velocity between earth and the instrument (=spacecraft).
  //
  //    a) Wavelengths: wl_spacecraft = wl_earth_rest_frame / (1.d0 + v_ei/c)
  //    b) Wavenumbers: wn_spacecraft = wn_earth_rest_frame * (1.d0 + v_ei/c)
  ArrayAd<double, 1> spec_dom_ad(Spec.spectral_domain().data_ad());
  if (spec_dom_ad.number_variable() == 0) {
    spec_dom_ad.resize_number_variable(coefficient().number_variable());
  }
  for(int w_idx = 0; w_idx < Spec.spectral_domain().data().rows(); w_idx++) {
    if (Spec.spectral_domain().type_preference() == SpectralDomain::PREFER_WAVELENGTH)
      spec_dom_ad(w_idx) = spec_dom_ad(w_idx) 
        / (1.0 + coefficient()(0) / OldConstant::speed_of_light.value);
    else
      spec_dom_ad(w_idx) = spec_dom_ad(w_idx) 
        * (1.0 + coefficient()(0) / OldConstant::speed_of_light.value);
  }
}

boost::shared_ptr<SpectrumEffect> InstrumentDoppler::clone() const
{
  return boost::shared_ptr<SpectrumEffect>(new InstrumentDoppler(DoubleWithUnit(coeff.value()(0), vel_units), used_flag(0)));
}

void InstrumentDoppler::print(std::ostream& Os) const
{
  Os << "InstrumentDoppler" << std::endl
     << "   Relative velocity: " << coefficient().value() << " " << vel_units << std::endl;
}
