// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "bad_sample_noise_model.h"
%}
%base_import(noise_model)

%fp_shared_ptr(FullPhysics::BadSampleNoiseModel);

namespace FullPhysics {
class BadSampleNoiseModel : public NoiseModel {
public:
  BadSampleNoiseModel(const boost::shared_ptr<NoiseModel>& 
		      Underlying_noise_model,
		      const blitz::Array<double, 2>& Bad_sample_mask,
                      double Bad_sample_uncer);
  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const;
  %python_attribute(bad_sample_uncertainty, double);
  %python_attribute(underlying_noise_model, 
		    const boost::shared_ptr<NoiseModel>&);
  %python_attribute(bad_sample_mask, 
		    const blitz::Array<double, 2>&);
};
}
