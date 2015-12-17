#ifndef NLLS_PROBLEM_STATE_H
#define NLLS_PROBLEM_STATE_H
#include <problem_state.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The state for a NLLS problem with implemented residual
///        function and its Jacobian
///
/// NLLSProblemState is used for the *Non-Linear Least Squares*
/// problem that its  residual function and Jacobian are implemented.
/// With this class one can store the current point in the parameter
/// space (the state), the value of the residual vector function at
/// that point, and the Jacobian matrix function at the same point.
///
/// Design related question(s):
///   -# Why is NLLSProblemState not derived from CostFuncDiffState?
///   -# Why not first implement a problem state class that only adds
///      residual to the state and then derive this class (with
///      residual and Jacobian) from that class (similar to
///      CostFuncState and CostFuncDiffState)?
///
/// Answers to the above question(s):
///   -# Keep in mind that the purpose of the classes in the class 
///      hierarchy rooted at ProblemState is to maintain computationally
///      expensive components of the cost function. If the residual 
///      and the Jacobian of a NLLS problem are computed, then the 
///      cost function and its gradient can be computed very fast.
///      If f(x) and J(x) are the residual and the Jacobian of the
///      NLLS problem respectively evaluated at x, then the cost 
///      function and its gradient respectively are
///      \f[
///          \frac{1}{2}\parallel f(x) \parallel^2
///      \f]
///      and
///      \f[
///          J(x)^T f(x)
///      \f]
///   -# If the Jacobian of a residual function is not available,
///      then the only methods that can solve the optimization
///      problem use the cost function.  In other words, in this 
///      case the problem is solved in a form that is presented by
///      CostFunc class.  For the problem in CostFunc form,
///      CostFuncState state class is sufficient.
//-----------------------------------------------------------------------

class NLLSProblemState : 
    virtual public ProblemState {

public:


//-----------------------------------------------------------------------
/// \brief Default constructor
//-----------------------------------------------------------------------

  NLLSProblemState() {}


//-----------------------------------------------------------------------
/// \brief Copy constructor
/// 
/// \param[in] s 
///            another NLLSProblemState
//-----------------------------------------------------------------------

  NLLSProblemState(const NLLSProblemState& s)
  { set(s); }


  virtual ~NLLSProblemState() {}


//-----------------------------------------------------------------------
/// \brief Makes self a copy of the input state
///
/// This method makes the object, for which it is called, a copy of
/// the input state.
/// 
/// \param[in] s 
///            another NLLSProblemState
//-----------------------------------------------------------------------

  virtual void set(const NLLSProblemState& s);


//-----------------------------------------------------------------------
/// \brief Deletes data contents.
///
/// This method deletes state.  If needed, it must be reimplemented
/// by other classes derived from this class to delete other saved 
/// components associated with the state as well.
//-----------------------------------------------------------------------

  virtual void clear()
  { ProblemState::clear(); R.free(); J.free(); }


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSProblemState"; }


protected:

  blitz::Array<double, 1> R;
  blitz::Array<double, 2> J;

};
}
#endif
