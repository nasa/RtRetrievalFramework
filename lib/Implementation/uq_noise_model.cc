#include "uq_noise_model.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(UqNoiseModel, NoiseModel)
.def(luabind::constructor<const HdfFile&,const HdfSoundingId&,const blitz::Array<double, 1>&>())
REGISTER_LUA_END()
#endif

void UqNoiseModel::read_hdf_noise(const HdfFile& Hfile, const HdfSoundingId& Sounding_id)
{
    Range ra = Range::all();

    TinyVector<int, 3> sz = Hfile.read_shape<3>("InstrumentHeader/snr_coef");
    Array<double, 3> snr_coefs = 
    Hfile.read_field<double, 3>("InstrumentHeader/snr_coef", 
       TinyVector<int, 3>(0,0,0), TinyVector<int, 3>(sz[0], sz[1], 2));

    coef_photon_.reference(snr_coefs(ra, ra, 0));
    coef_background_.reference(snr_coefs(ra, ra, 1));

    // Reorder dimensions as num_pixels x num_spectrometers
    coef_photon_.transposeSelf(secondDim, firstDim);
    coef_background_.transposeSelf(secondDim, firstDim);
}


