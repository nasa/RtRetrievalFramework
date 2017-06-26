#include "spectrally_resolved_noise.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(SpectrallyResolvedNoise, NoiseModel)
.def(luabind::constructor<boost::shared_ptr<NoiseModel> &>())
.def("set_noise_coefficients", &SpectrallyResolvedNoise::set_noise_coefficients)
REGISTER_LUA_END()
#endif
  
//-----------------------------------------------------------------------
/// Add noise coefficents for a given band index. The coefficeints
/// per band do not need to be added in any particular order
//-----------------------------------------------------------------------

void SpectrallyResolvedNoise::set_noise_coefficients(int Spec_index, const blitz::Array<double, 1> Noise_coeff) {
  if(band_noise_coeffs_.size() <= Spec_index)
    band_noise_coeffs_.resize(Spec_index+1);
  band_noise_coeffs_[Spec_index].reference(Noise_coeff);
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
  uncert(all) = uncert(all) * band_noise_coeffs_[Spec_index](all);
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
