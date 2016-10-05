// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_max_a_posteriori.h"
%}
%base_import(nlls_problem)
%base_import(nlls_problem_state)
%import "max_a_posteriori.i"
%fp_shared_ptr(FullPhysics::NLLSMaxAPosteriori);

namespace FullPhysics {
class NLLSMaxAPosteriori: public NLLSProblem, public NLLSProblemState {
public:
  NLLSMaxAPosteriori(const boost::shared_ptr<MaxAPosteriori>& map, bool together=false);
  virtual ~NLLSMaxAPosteriori();
  %python_attribute_nonconst(residual, blitz::Array<double, 1>)
  %python_attribute_nonconst(jacobian, blitz::Array<double, 2>)
  %python_attribute(residual_size, int)
  %python_attribute(expected_parameter_size, int)
  %python_attribute_with_set(parameters, blitz::Array<double, 1>)
  %python_attribute_nonconst(max_a_posteriori, boost::shared_ptr<MaxAPosteriori>)
};
}
