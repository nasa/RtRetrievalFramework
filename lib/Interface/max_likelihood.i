// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "max_likelihood.h"
%}
%base_import(model_measure)
%fp_shared_ptr(FullPhysics::MaxLikelihood);

namespace FullPhysics {
class MaxLikelihood : virtual public ModelMeasure {
public:
  MaxLikelihood();
  virtual ~MaxLikelihood();
};
}
