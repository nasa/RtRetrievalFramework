// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "oco_noise_model.h"
#include "oco_sounding_id.h"
#include "hdf_file.h"
#include "state_vector.h"
#include "instrument.h"
%}
%base_import(noise_model)
%import "hdf_file.i"
%import "oco_sounding_id.i"
%import "instrument.i"
%fp_shared_ptr(FullPhysics::OcoNoiseModel);

namespace FullPhysics {
class OcoNoiseModel: public NoiseModel {
public:
  OcoNoiseModel(const HdfFile& Hfile, const HdfSoundingId& Sounding_id);
  OcoNoiseModel(const HdfFile& Hfile, const HdfSoundingId& Sounding_id, const blitz::Array<double, 1>& Max_meas_signal);
  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const;
};
}
