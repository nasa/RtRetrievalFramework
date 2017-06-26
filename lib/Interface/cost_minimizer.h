#ifndef COST_MINIMIZER_H
#define COST_MINIMIZER_H
#include <iterative_solver.h>
#include <cost_func.h>


namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for all iterative cost minimizers
///        that do not require derivatives of any order.
///
/// This is the base class for methods that iteratively
/// minimize a scalar cost function.  In other words, this 
/// is a base class for methods that find a point in the
/// parameter space where the cost function is at least 
/// locally minimum.
///
/// This class is associated with a problem (CostFunc)
/// because the problem interface is determined:
///   - provide a point in the parameter space
///   - evaluate the cost function (a scalar function)
///     at the point
//-----------------------------------------------------------------------

class CostMinimizer : 
    public IterativeSolver {

public:


//-----------------------------------------------------------------------
/// \brief Constructor
/// 
/// \param[in] max_cost_function_calls 
///            read related base class comments
///
/// \param[in] dx_tol_abs
///            read related base class comments
///
/// \param[in] dx_tol_rel
///            read related base class comments
///
/// \param[in] p
///            The cost minimization problem
///
/// \param[in] vrbs
///            read related base class comments
//-----------------------------------------------------------------------

  CostMinimizer(int max_cost_function_calls, 
                double dx_tol_abs, double dx_tol_rel, 
                const boost::shared_ptr<CostFunc>& p,
                bool vrbs)
    : IterativeSolver(max_cost_function_calls, dx_tol_abs, dx_tol_rel, vrbs),
      P(p)
  {}


  virtual ~CostMinimizer() {}


//-----------------------------------------------------------------------
/// \brief Returns the cost minimization problem
///
/// This method returns the cost minimization problem that
/// is passed to the constructor of the solver.
///
/// \return Cost-function problem
//-----------------------------------------------------------------------

  const boost::shared_ptr<CostFunc>& cost_min_problem() const
  { return P; }


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "CostMinimizer"; }


protected:

  boost::shared_ptr<CostFunc> P;

};
}
#endif
