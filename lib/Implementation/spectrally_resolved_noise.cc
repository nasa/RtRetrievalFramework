#include "spectrally_resolved_noise.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(SpectrallyResolvedNoise, NoiseModel)
.def(luabind::constructor<boost::shared_ptr<NoiseModel> &>())
.def("set_full_noise_scaling", &SpectrallyResolvedNoise::set_full_noise_scaling)
.def("set_single_noise_scaling", &SpectrallyResolvedNoise::set_single_noise_scaling)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create a new SpectrallyResolvedNoise modifying an existing noise
/// model. The created object should be used in place of the base model.
/// The noise coefficents for each band must be set individually since
/// they might vary in size
//-----------------------------------------------------------------------
SpectrallyResolvedNoise::SpectrallyResolvedNoise(boost::shared_ptr<NoiseModel> Base_model) 
: base_model_(Base_model)
{
}
  
//-----------------------------------------------------------------------
/// Add noise scalings for a given band index where each array value
/// corresponds to a sample index in the base noise uncertainty.
//-----------------------------------------------------------------------

void SpectrallyResolvedNoise::set_full_noise_scaling(int Spec_index, const blitz::Array<double, 1> Noise_scaling) {
  if(band_noise_coeffs_.size() <= Spec_index) {
    band_noise_coeffs_.resize(Spec_index+1);
    band_single_value_.resize(Spec_index+1);
  }
  band_noise_coeffs_[Spec_index].reference(Noise_scaling);
  band_single_value_[Spec_index] = false;
}

//-----------------------------------------------------------------------
/// Add noise scalings for a given band index where a single value applies
/// across the whole band.
//-----------------------------------------------------------------------

void SpectrallyResolvedNoise::set_single_noise_scaling(int Spec_index, double Noise_scaling) {
  if(band_noise_coeffs_.size() <= Spec_index) {
    band_noise_coeffs_.resize(Spec_index+1);
    band_single_value_.resize(Spec_index+1);
  }
  band_noise_coeffs_[Spec_index].resize(1);
  band_noise_coeffs_[Spec_index](0) = Noise_scaling;
  band_single_value_[Spec_index] = true;
}

//-----------------------------------------------------------------------
/// Applies spectrally resolved noise coefficents to base model's
/// uncertainty. If coefficents are not set for a a spec_index then
/// the base model's uncertainty is returned as is.
//-----------------------------------------------------------------------

blitz::Array<double, 1> SpectrallyResolvedNoise::uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const {
  Array<double, 1> uncert = base_model_->uncertainty(Spec_index, Radiance);
  if (band_noise_coeffs_.size() <= Spec_index or band_noise_coeffs_[Spec_index].rows() == 0)
    return uncert;

  Range all = Range::all(); 
  if(band_single_value_[Spec_index]) {
    uncert(all) = uncert(all) * band_noise_coeffs_[Spec_index](0);
  } else {
    uncert(all) = uncert(all) * band_noise_coeffs_[Spec_index](all);
  }
  return uncert;
}

void SpectrallyResolvedNoise::print(std::ostream& Os) const
{
  Os << "Spectrally Resolved Noise Modification:" << std::endl;
  OstreamPad opad(Os, "    ");
  opad << "Coefficient sizes:" << std::endl;
  for(int spec_idx = 0; spec_idx < band_noise_coeffs_.size(); spec_idx++)
    opad << "  [" << spec_idx << "]: " << band_noise_coeffs_[spec_idx].rows() << std::endl;
  opad.strict_sync();
}
