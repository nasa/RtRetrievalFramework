#ifndef FM_NLLS_PROBLEM
#define FM_NLLS_PROBLEM
#include <nlls_problem.h>
#include "state_vector.h"
#include "forward_model.h"
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/****************************************************************//**
  This is the forward model cost function, rewritten to match the
  standard format used in nonlinear least squares solvers.
*******************************************************************/
class FmNLLSProblem : public NLLSProblem {
public:
  FmNLLSProblem(const boost::shared_ptr<ForwardModel>& Fm,
   const blitz::Array<double, 1>& Rad,
   const blitz::Array<double, 1>& Rad_uncer,
   const blitz::Array<double, 1> X_apriori,
   const blitz::Array<double, 2> Apriori_cov);
  virtual ~FmNLLSProblem() {}
  virtual int residual_size() const { return rad.rows() + x_a.rows(); }
  virtual int expected_parameter_size() const { return x_a.rows(); }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "FmNLLSProblem"; }
private:
  boost::shared_ptr<ForwardModel> fm;
  blitz::Array<double, 1> rad;
  blitz::Array<double, 1> x_a;
  blitz::Array<double, 1> se_sqrt_inv;
  blitz::Array<double, 2> sa_sqrt_inv;
};
}

#endif
