#ifndef UQ_NOISE_MODEL_H
#define UQ_NOISE_MODEL_H
#include "oco_noise_model.h"

namespace FullPhysics {
/****************************************************************//**
  This class creates a UQ noise model using inputs from the supplied file
*******************************************************************/
class UqNoiseModel: public OcoNoiseModel {
public:

    /// Creates a new UQ noise model from the L1B HDF5 file
    UqNoiseModel(const HdfFile& Hfile, const HdfSoundingId& Sounding_id, const blitz::Array<double, 1>& Max_meas_signal) {
        max_ms_.resize(Max_meas_signal.shape());
        max_ms_ = Max_meas_signal;

        read_hdf_noise(Hfile, Sounding_id);
    }

protected:
    virtual void read_hdf_noise(const HdfFile& Hfile, const HdfSoundingId& Sounding_id);
};
}
#endif


