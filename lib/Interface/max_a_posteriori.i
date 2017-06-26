// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "max_a_posteriori.h"
%}
%base_import(model_measure)
%fp_shared_ptr(FullPhysics::MaxAPosteriori);

namespace FullPhysics {
class MaxAPosteriori : virtual public ModelMeasure {
public:
  MaxAPosteriori(const blitz::Array<double, 1>& a_priori_params,
                 const blitz::Array<double, 2>& a_priori_cov);
  virtual ~MaxAPosteriori();
  %python_attribute(a_priori_params, blitz::Array<double, 1>)
  %python_attribute(a_priori_cov, blitz::Array<double, 2>);

  %python_attribute(param_a_priori_uncertainty, blitz::Array<double, 1>)
  %python_attribute(parameter_a_priori_diff, blitz::Array<double, 1>)
  %python_attribute(cov_weighted_parameter_a_priori_diff, blitz::Array<double, 1>)
  %python_attribute(a_priori_cov_chol_inv, blitz::Array<double, 2>)
  %python_attribute_nonconst(weighted_model_measure_diff_aug, 
			     blitz::Array<double, 1>)
  %python_attribute_nonconst(a_posteriori_covariance, blitz::Array<double, 2>)
  %python_attribute(a_priori_cov_chol, blitz::Array<double, 2>)
  %python_attribute_nonconst(param_a_posteriori_uncertainty, blitz::Array<double, 1>)
};
}
