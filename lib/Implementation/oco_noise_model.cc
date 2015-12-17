#include "oco_noise_model.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(OcoNoiseModel, NoiseModel)
.def(luabind::constructor<const HdfFile&,const HdfSoundingId&,const blitz::Array<double, 1>&>())
REGISTER_LUA_END()
#endif

void OcoNoiseModel::read_hdf_noise(const HdfFile& Hfile, 
				   const HdfSoundingId& Sounding_id) 
{
  int sindex = Sounding_id.sounding_number(); 
  Range ra = Range::all();

  TinyVector<int, 4> sz = Hfile.read_shape<4>("InstrumentHeader/snr_coef");
  Array<double, 4> snr_coefs = 
    Hfile.read_field<double, 4>("InstrumentHeader/snr_coef", 
				TinyVector<int, 4>(0,sindex,0,0),
				TinyVector<int, 4>(sz[0], 1, sz[2], 2));
  coef_photon_.reference(snr_coefs(ra, 0, ra, 0));
  coef_background_.reference(snr_coefs(ra, 0, ra, 1));

  // Reorder dimensions as num_pixels x num_spectrometers
  coef_photon_.transposeSelf(secondDim, firstDim);
  coef_background_.transposeSelf(secondDim, firstDim);
}

//-----------------------------------------------------------------------
/// Calculate uncertainty for OCO using SNR coefficient tables.
//-----------------------------------------------------------------------

blitz::Array<double, 1> OcoNoiseModel::uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const {
  range_max_check(Spec_index, coef_photon_.extent(secondDim)+1);
  
  if (Radiance.extent(firstDim) != coef_photon_.extent(firstDim)) {
    Exception err;
    err << "Size of radiance array: " << Radiance.extent(firstDim) << " "
	<< "does not match size of SNR coefs data for the pixel dimension: " 
	<< coef_photon_.extent(firstDim);
    throw err;
  }
  
  Range ra = Range::all();
  Array<double,1> res(Radiance.extent(firstDim));
  res = (100.0 * blitz::where(Radiance > 0, Radiance, 0.0) / max_ms_(Spec_index)) * coef_photon_(ra,Spec_index) * coef_photon_(ra,Spec_index);
  res = res + coef_background_(ra,Spec_index) * coef_background_(ra,Spec_index);
  res = sqrt(res);
  res *= max_ms_(Spec_index) / 100.0;

  return res;
}

void OcoNoiseModel::print(std::ostream& Os) const
{
  Os << "Oco Noise Model:" << std::endl
     << "  Num Bands:    " << coef_photon_.cols() << std::endl
     << "  Num Colors:   " << coef_photon_.rows() << std::endl;
}
