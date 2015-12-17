// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_max_likelihood.h"
%}
%base_import(nlls_problem)
%base_import(nlls_problem_state)
%import "max_likelihood.i"
%fp_shared_ptr(FullPhysics::NLLSMaxLikelihood);

namespace FullPhysics {
class NLLSMaxLikelihood : public NLLSProblem, public NLLSProblemState {
public:
  NLLSMaxLikelihood(const boost::shared_ptr<MaxLikelihood>& ml, bool together=false);
  virtual ~NLLSMaxLikelihood();
  %python_attribute_nonconst(residual, blitz::Array<double, 1>)
  %python_attribute_nonconst(jacobian, blitz::Array<double, 2>)
  %python_attribute(residual_size, int)
  %python_attribute(expected_parameter_size, int)
  %python_attribute_with_set(parameters, blitz::Array<double, 1>)
  %python_attribute_nonconst(max_likelihood, boost::shared_ptr<MaxLikelihood>)
};
}
