#include "precomputed_noise_model.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(PrecomputedNoiseModel, NoiseModel)
.def(luabind::constructor<const HeritageFile&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Reads noise model data directly from a file and passes through as is.
//-----------------------------------------------------------------------

PrecomputedNoiseModel::PrecomputedNoiseModel(const HeritageFile& Noise_file) {
  
  // Read noise column, failing if not present
  Array<double, 1> noise_column = Noise_file.data("Noise");

  // Read indexing of pixels
  std::vector<int> start_pixels = Noise_file.value<std::vector<int> >("start_pixels");

  file_noise_values.resize(start_pixels.size());

  for(unsigned int spec_idx = 0; spec_idx < start_pixels.size(); spec_idx++) {
    Range read_range;
    // -1 for fortran -> c indexing
    if(spec_idx != start_pixels.size()-1) {
      read_range = Range(start_pixels[spec_idx]-1, start_pixels[spec_idx+1]-2);
    } else {
      read_range = Range(start_pixels[spec_idx]-1, noise_column.rows()-1);
    }
    file_noise_values[spec_idx].resize(read_range.length());
    file_noise_values[spec_idx](Range::all()) = noise_column(read_range);
  }
  
}

//-----------------------------------------------------------------------
/// Returns uncertainty values read from a file. In this instance the
/// radiance values are not used for the purpose of the noise computation
//-----------------------------------------------------------------------

blitz::Array<double, 1> PrecomputedNoiseModel::uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const {
  Array<double, 1> result = file_noise_values[Spec_index];

  if (result.rows() != Radiance.rows()) {
    std::stringstream err_msg;
    err_msg << "Size of noise data: "
	    << result.rows()
	    << " for spectrometer: " << Spec_index
	    << " does not match size expected from radiance data: "
	    << Radiance.rows();
    throw Exception(err_msg.str());
  }

  return result;
}

void PrecomputedNoiseModel::print(std::ostream& Os) const
{
  Os << "Precomputed Noise Model:" << std::endl
     << "  Num Bands:    " << file_noise_values.size() << std::endl;
}

