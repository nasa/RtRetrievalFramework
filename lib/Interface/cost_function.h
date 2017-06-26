#ifndef COST_FUNCTION_H
#define COST_FUNCTION_H
#include "printable.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This class calculates a cost function, along with a jacobian.
*******************************************************************/
class CostFunction : public Printable<CostFunction> {
public:
  virtual ~CostFunction() {}

//-----------------------------------------------------------------------
/// For the given value of X, calculate the residuals and jacobians.
///
/// The residual is defined as F(x) - y, so the Jacobian of the
/// residuals is the same as the Jacobian of F(x).
///
/// \param X Input value
/// \param Residual On exit, set to the residual of the cost function
/// \param Se On exit, the covariance of the residual. We assume the 
///     covariance is a diagonal matrix, and just return the diagonal
///     elements
/// \param Jacobian On exit, the Jacobian of the cost function.
//-----------------------------------------------------------------------

  virtual void cost_function(const blitz::Array<double, 1>& X,
			blitz::Array<double, 1>& Residual,
			blitz::Array<double, 1>& Se,
			blitz::Array<double, 2>& Jacobian) const = 0;
			
//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const {Os << "CostFunction";}
};
}
#endif
