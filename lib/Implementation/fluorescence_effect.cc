#include "fluorescence_effect.h"
#include "forward_model_spectral_grid.h"
#include "ostream_pad.h"

#include <boost/progress.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(FluorescenceEffect, SpectrumEffect)
.def(luabind::constructor<const blitz::Array<double, 1>&,
                          const blitz::Array<bool, 1>&,
                          const boost::shared_ptr<RtAtmosphere>&,
                          const boost::shared_ptr<StokesCoefficient>&,
                          const DoubleWithUnit&, 
                          const int,
                          const DoubleWithUnit&,
                          const Unit&>())
REGISTER_LUA_END()
#endif
  
FluorescenceEffect::FluorescenceEffect
(const blitz::Array<double, 1>& Coeff,
 const blitz::Array<bool, 1>& Used_flag,
 const boost::shared_ptr<RtAtmosphere>& Atm,
 const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
 const DoubleWithUnit& Lza, 
 const int Spec_index,
 const DoubleWithUnit& Reference,
 const Unit& Retrieval_unit) 
: SpectrumEffectImpBase(Coeff, Used_flag),
  lza(Lza), reference(Reference), retrieval_unit(Retrieval_unit), 
  spec_index(Spec_index), stokes_coef(Stokes_coef)
{
  // We need to use the AtmosphereOco specific interfaces
  atm_oco = boost::dynamic_pointer_cast<AtmosphereOco>(Atm);

  // Use a map cache, because due to LSI and non-uniform sampling we can
  // not count on all values being cached that we need in the right order
  atm_oco->column_optical_depth_cache().reset( new ArrayAdMapCache<double, double, 1>() );
}

void FluorescenceEffect::apply_effect(Spectrum& Spec,
	      const ForwardModelSpectralGrid& Forward_model_grid) const
{
  AccumulatedTimer tm("FluorescenceEffect");
  {
    FunctionTimer ft = tm.function_timer();

    double lza_sec = 1.0 / cos(lza.convert(units::rad).value);
    double reference_wn = reference.convert_wave(units::inv_cm).value;

    Logger::info() << "Adding fluorescence to band: " + boost::lexical_cast<std::string>(spec_index + 1) + "\n";

    // Radiance array of spectrum to add fluoresence to
    ArrayAd<double, 1>& rt_solar_rad = Spec.spectral_range().data_ad();

    // Access fluoresence retireved values in the correct units
    double conv_factor = FullPhysics::conversion(retrieval_unit, Spec.spectral_range().units());
    AutoDerivative<double> fs_ref = coefficient()(0) * conv_factor;
    AutoDerivative<double> slope = coefficient()(1);

    // Grid to calculate on and Spectrum class to save into for use in interpolation
    // Must use same grid as used in RT and interpolated later (if needed) so that points
    // match up correctly when doing final addition
    Array<double, 1> spec_samp_wn =
      Forward_model_grid.high_resolution_grid(spec_index).wavenumber();
    ArrayAd<double, 1> fluor_effect_rt;

    // Calculate fluorescence contribution that matches what would have been calculated if done inside RT
    for(int wn_idx = 0; wn_idx < spec_samp_wn.rows(); wn_idx++) {
      AutoDerivative<double> o2_col_abs = atm_oco->column_optical_depth(spec_samp_wn(wn_idx), spec_index, "O2");

      AutoDerivative<double> f_surf = fs_ref * (1.0 + slope * (spec_samp_wn(wn_idx) - reference_wn));
      if(wn_idx ==0) {
          fluor_effect_rt.resize(spec_samp_wn.rows(),
                                 std::max(f_surf.number_variable(), o2_col_abs.number_variable()));
      }
      fluor_effect_rt(wn_idx) = f_surf * exp(-1 * o2_col_abs * lza_sec);     
    }
    
    Spectrum f_contrib_spectrum(Forward_model_grid.high_resolution_grid(spec_index), SpectralRange(fluor_effect_rt, Spec.spectral_range().units()));
    
    // Interpolate to obtain contribution we should add to rt radiances which have been interpolated and solar
    // doppler shifted
    Spectrum fluor_effect_interpolated = Forward_model_grid.interpolate_spectrum(f_contrib_spectrum, spec_index);

    // Save interpolated fluoresence contribution for use in output products
    f_contrib_ad.reference(fluor_effect_interpolated.spectral_range().data_ad());

    // Check that interpolated fluorescence and radiances agree in size
    if (rt_solar_rad.rows() != f_contrib_ad.rows()) {
      std::stringstream err_msg;
      err_msg << "Interpolated fluorescence size: " << f_contrib_ad.rows() 
              << " does not match radiance size: " << rt_solar_rad.rows();
      throw Exception(err_msg.str());
    }

    // Add fluorescence contribution to radiances
    if(rt_solar_rad.number_variable() == 0 && f_contrib_ad.number_variable() > 0) {
      rt_solar_rad.resize_number_variable(f_contrib_ad.number_variable());
    }

    for(int wn_idx = 0; wn_idx < rt_solar_rad.rows(); wn_idx++) {
      rt_solar_rad(wn_idx) = rt_solar_rad(wn_idx) + 
	stokes_coef->stokes_coefficient()(spec_index, 0) * f_contrib_ad(wn_idx);
    }
  }
  Logger::info() << tm << "\n";
}

boost::shared_ptr<SpectrumEffect> FluorescenceEffect::clone() const
{
  return boost::shared_ptr<SpectrumEffect>(new FluorescenceEffect(coeff.value(), used_flag,
                                           atm_oco, stokes_coef, lza, spec_index,
                                           reference, retrieval_unit));
}

void FluorescenceEffect::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "FluorescenceEffect" << std::endl
     << "  Coefficient:     " << coefficient().value() << "\n"
     << "  Retrieval flag:  " << used_flag_value() << "\n"
     << "  Zenith angle:    " << lza << "\n"
     << "  Reference point: " << reference << "\n";
  opad << *stokes_coef;
  opad.strict_sync();
}
