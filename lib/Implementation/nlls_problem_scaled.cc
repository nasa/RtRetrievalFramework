#include "nlls_problem_scaled.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<CostFunc> create_nlls_problem_scaled
(const Array<double, 1>& S,
 const boost::shared_ptr<CostFunc>& Cf)
{
  boost::shared_ptr<NLLSProblem> p = boost::dynamic_pointer_cast<NLLSProblem>(Cf);
  return boost::shared_ptr<CostFunc>(new NLLSProblemScaled(S, p));
}
			   
REGISTER_LUA_DERIVED_CLASS(NLLSProblemScaled, CostFunc)
.def(luabind::constructor< const Array<double, 1>&, const boost::shared_ptr<NLLSProblem>& >())
.def("parameters", (void( NLLSProblemScaled::*)(const blitz::Array<double, 1>&))&NLLSProblemScaled::parameters)
.def("scale_parameters", &NLLSProblemScaled::scale_parameters)
.scope
[
 luabind::def("create", &create_nlls_problem_scaled)
]
REGISTER_LUA_END()
#endif



NLLSProblemScaled::NLLSProblemScaled( const Array<double, 1>& s,
                   const boost::shared_ptr<NLLSProblem>& p ) 
  : NLLSProblem(), S(s.copy()), P(p)
{
  if(s.rows() != P->expected_parameter_size()) {
    Exception e;
    e << "Parameter space scaling diagonal matrix and parameters not equal in size:\n"
      << " Parameter size " << P->expected_parameter_size()
      << " Parameter space scaling diagonal matrix size " << s.rows();
    throw e;
  }
}


Array<double, 1> NLLSProblemScaled::residual()
{
  assert_parameter_set_correctly();
  P->parameters(unscale_parameters(X));
  Array<double, 1> R(P->residual());

  //  Depending on the developers intention, one of
  //  the following two methods can be used to update
  //  the cost evaluation counter.  Just be consistent
  //  with the method used to update gradient (first 
  //  order derivatives) evaluation counter in the
  //  method NLLSProblemScaled::jacobian()
  //
  set_num_cost_evaluations(P->num_residual_evaluations());
//  increment_num_cost_evaluations();

  return R;
}


Array<double, 2> NLLSProblemScaled::jacobian()
{
  assert_parameter_set_correctly();
  P->parameters(unscale_parameters(X));
  Array<double, 2> J(P->jacobian());
  firstIndex i1; secondIndex i2;
  J = J(i1,i2)*S(i2);

  //  Depending on the developers intention, one of
  //  the following two methods can be used to update
  //  the gradient (first order derivatives) evaluation
  //  counter.  Just be consistent with the method used
  //  to update cost evaluation counter in the method
  //  NLLSProblemScaled::residual()
  //
  set_num_der1_evaluations(P->num_jacobian_evaluations());
//  increment_num_der1_evaluations();

  return J;
}


Array<double, 1> NLLSProblemScaled::scale_parameters(const Array<double, 1>& x) const
{
  Array<double,1> _x_(x.copy());
  for(int i=0; i<_x_.rows(); i++)
    if(S(i)) _x_(i) /= S(i);
  return _x_;
}


Array<double, 1> NLLSProblemScaled::unscale_parameters(const Array<double, 1>& x) const
{
  Array<double,1> _x_(x.copy());
  for(int i=0; i<_x_.rows(); i++)
    if(S(i)) _x_(i) *= S(i);
  return _x_;
}


void NLLSProblemScaled::parameters(const blitz::Array<double, 1>& x)
{
  P->parameters(unscale_parameters(x));
  NLLSProblem::parameters(x);
}
