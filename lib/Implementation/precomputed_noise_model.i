// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "precomputed_noise_model.h"
%}
%base_import(noise_model)
%import "heritage_file.i"

%fp_shared_ptr(FullPhysics::PrecomputedNoiseModel);

namespace FullPhysics {
class PrecomputedNoiseModel: public NoiseModel {
public:
  PrecomputedNoiseModel(const HeritageFile& Noise_file);
  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const;
};
}
