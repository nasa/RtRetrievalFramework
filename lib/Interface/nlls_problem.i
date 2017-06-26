// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_problem.h"
%}
%base_import(cost_func_diff)

%fp_shared_ptr(FullPhysics::NLLSProblem);

namespace FullPhysics {
class NLLSProblem : public CostFuncDiff {
public:
  NLLSProblem();
  virtual ~NLLSProblem();
  %python_attribute_nonconst(cost, double);
  %python_attribute_nonconst(residual, blitz::Array<double, 1>);
  virtual blitz::Array<double, 1> residual_x(const blitz::Array<double, 1>& x);
  %python_attribute_nonconst(jacobian, blitz::Array<double, 2>);
  virtual blitz::Array<double, 2> jacobian_x(const blitz::Array<double, 1>& x);
  virtual void residual_jacobian(
    blitz::Array<double, 1>& OUTPUT, blitz::Array<double, 2>& OUTPUT);
  virtual void residual_jacobian_x(const blitz::Array<double, 1>& x,
    blitz::Array<double, 1>& OUTPUT, blitz::Array<double, 2>& OUTPUT);
  %python_attribute(num_residual_evaluations, int);
  %python_attribute(num_jacobian_evaluations, int);
  %python_attribute_abstract(residual_size, int);
};
}
