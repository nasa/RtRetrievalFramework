// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_problem_scaled.h"
%}
%base_import(nlls_problem)
%fp_shared_ptr(FullPhysics::NLLSProblemScaled);

namespace FullPhysics {
class NLLSProblemScaled : public NLLSProblem {
public:
  NLLSProblemScaled( const blitz::Array<double, 1>& s,
                     const boost::shared_ptr<NLLSProblem>& p );
  virtual ~NLLSProblemScaled();
  %python_attribute_nonconst(residual, blitz::Array<double, 1>)
  %python_attribute_nonconst(jacobian, blitz::Array<double, 2>)
  %python_attribute(residual_size, int)
  %python_attribute(expected_parameter_size, int)
  %python_attribute_with_set(parameters, blitz::Array<double, 1>)
  %python_attribute_nonconst(nlls_problem, boost::shared_ptr<NLLSProblem>)
  virtual blitz::Array<double, 1> scale_parameters(const blitz::Array<double, 1>& x) const;
  virtual blitz::Array<double, 1> unscale_parameters(const blitz::Array<double, 1>& x) const;
};
}
