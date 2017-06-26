// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "max_likelihood_oco.h"
%}
%base_import(max_likelihood)
%base_import(model_measure_oco)
%import "forward_model.i"
%fp_shared_ptr(FullPhysics::MaxLikelihoodOCO);

namespace FullPhysics {
class MaxLikelihoodOCO : public MaxLikelihood, public ModelMeasureOCO {
public:
  MaxLikelihoodOCO(const boost::shared_ptr<ForwardModel>& fm);
  virtual ~MaxLikelihoodOCO();
};
}
