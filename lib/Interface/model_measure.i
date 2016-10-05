// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "model_measure.h"
%}
%base_import(model_state)
%fp_shared_ptr(FullPhysics::ModelMeasure);

namespace FullPhysics {
class ModelMeasure : public ModelState {
public:
  ModelMeasure(const blitz::Array<double, 1>& measurement, 
               const blitz::Array<double, 1>& measurement_error_cov);
  ModelMeasure();
  virtual ~ModelMeasure();
  virtual void model_eval() = 0;
  %python_attribute_nonconst(model, blitz::Array<double, 1>)
  virtual blitz::Array<double, 1> model_x(const blitz::Array<double, 1>& x);
  virtual void jacobian_eval() = 0;
  %python_attribute_nonconst(jacobian, blitz::Array<double, 2>)
  virtual blitz::Array<double, 2> jacobian_x(const blitz::Array<double, 1>& x);
  virtual void model_jacobian_eval();
  virtual void model_jacobian(blitz::Array<double, 1>& m, 
			      blitz::Array<double, 2>& k);
  virtual void model_jacobian_x(const blitz::Array<double, 1>& x,
     blitz::Array<double, 1>& m, blitz::Array<double, 2>& k);
  %python_attribute(measurement, blitz::Array<double, 1>)
  %python_attribute(measurement_error_cov, blitz::Array<double, 1>)
  %python_attribute(measurement_size, int)
  virtual void assert_model_correct(const blitz::Array<double, 1>& m) const;
  virtual void assert_jacobian_correct(const blitz::Array<double, 2>& k) const;

  %python_attribute_nonconst(model_measure_diff, blitz::Array<double, 1>)
  %python_attribute_nonconst(uncert_weighted_model_measure_diff, 
			     blitz::Array<double, 1>)
  %python_attribute_nonconst(uncert_weighted_jacobian,
			     blitz::Array<double, 2>)
  %python_attribute_nonconst(uncert_weighted_jac_inner_product,
			     blitz::Array<double, 2>)
};
}
