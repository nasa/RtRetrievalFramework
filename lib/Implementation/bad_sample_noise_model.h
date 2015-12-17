#ifndef BAD_SAMPLE_NOISE_MODEL_H
#define BAD_SAMPLE_NOISE_MODEL_H
#include "noise_model.h"
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/****************************************************************//**
  When we have bad samples, we usually pass this to the spectral
  window to prevent the sample from even being used. However, it can
  be useful to instead include the bad sample but give it a really 
  large noise value (e.g., you want the residual calculated for a ARP
  investigation, but don't want it to have any weight in the
  retrieval). This adapter class takes a underlying noise model, but
  changes the uncertainty for any value that is marked as a bad sample.
*******************************************************************/
class BadSampleNoiseModel: public NoiseModel {
public:
  template <class T> 
  BadSampleNoiseModel(const boost::shared_ptr<NoiseModel>& 
		      Underlying_noise_model,
		      const blitz::Array<T, 2>& Bad_sample_mask,
                      double Bad_sample_uncer)
    : underlying_noise_model_(Underlying_noise_model),
      bad_sample_uncer_(Bad_sample_uncer)
  {
    using namespace blitz;
    bad_sample_mask_.resize(Bad_sample_mask.shape());
    for(int i = 0; i < bad_sample_mask_.rows(); ++i)
      bad_sample_mask_(i, Range::all()) = where(Bad_sample_mask(i, Range::all()), true, false);
  }
  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const;

  virtual void print(std::ostream& Os) const;

  double bad_sample_uncertainty() const { return bad_sample_uncer_; }
  const boost::shared_ptr<NoiseModel>& underlying_noise_model() const
  { return underlying_noise_model_;}
  const blitz::Array<bool, 2>& bad_sample_mask() const 
  {return bad_sample_mask_;}
private:
  boost::shared_ptr<NoiseModel> underlying_noise_model_;
  blitz::Array<bool, 2> bad_sample_mask_;
  double bad_sample_uncer_;
};
}
#endif


