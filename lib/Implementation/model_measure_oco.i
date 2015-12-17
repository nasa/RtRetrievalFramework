// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "model_measure_oco.h"
%}
%base_import(model_measure)
%fp_shared_ptr(FullPhysics::ModelMeasureOCO);

namespace FullPhysics {
class ModelMeasureOCO : virtual public ModelMeasure {
public:
  virtual ~ModelMeasureOCO();
  virtual void model_eval();
  virtual void jacobian_eval();
  virtual void model_jacobian_eval();
  %python_attribute(expected_parameter_size, int)
  %python_attribute_with_set(parameters,blitz::Array<double, 1>)
protected:
  ModelMeasureOCO();
};
}
