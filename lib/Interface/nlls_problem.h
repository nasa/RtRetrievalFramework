#ifndef NLLS_PROBLEM_H
#define NLLS_PROBLEM_H
#include <cost_func_diff.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for the Non-Linear Least Squares problem.
///
/// The class NLLSProblem is the base class for all problem classes
/// that implement a Non-Linear Least Squares (NLLS) problem.  The
/// two main interface components provided by this class are
///   - the residual of the problem (a vector function)
///   - the Jacobian of the residual (the first order derivatives
///     and a matrix function)
///
/// A NLLS problem can be solved by NLLS solvers such as the 
/// Gauss-Newton or the Levenberg-Marquardt methods.
///
/// NLLSProblem implements cost and gradient from its inherited
/// classes; therefore, an NLLS problem can also be solved by methods
/// that solve optimization problems of CostFunc or CostFuncDiff
/// form.
///
/// A DESIGN RELATED QUESTION:
///
/// Similar to CostFunc and CostFuncDiff, why don't we have a 
/// class that only adds the residual to the class interface
/// and then derive this class from that one to add the Jacobian
/// to the interface as well?
///
/// ANSWER TO THE ABOVE QUESTION:
///
/// As mentioned in the comments of CostFunc, CostFuncDiff and
/// NLLSProblem (this) classes, the problems of the forms
///   - cost only (CostFunc)
///   - cost and gradient (CostFuncDiff)
///   - residual and Jacobian (NLLSProblem)
///
/// can be solved by certain algorithms.  However, there is
/// no method for solving a problem of the form 
///   - residual only
///
/// If the first order derivatives (the Jacobian) of the residual
/// are not available, then the problem can be solved only by the
/// methods that solve a problem of cost-only (CostFunc) form 
/// after converting the residual vector function into a cost 
/// scalar function using the equation
/// \f[
///     c(x) = \frac{1}{2}\parallel f(x) \parallel^2
/// \f]
//-----------------------------------------------------------------------

class NLLSProblem : 
    public CostFuncDiff {

public:


//-----------------------------------------------------------------------
/// \brief Default Constructor
//-----------------------------------------------------------------------

  NLLSProblem() 
    : CostFuncDiff()
  {}


  virtual ~NLLSProblem() {}


//-----------------------------------------------------------------------
/// \brief Read comments on CostFunc::cost()
//-----------------------------------------------------------------------

  virtual double cost();


//-----------------------------------------------------------------------
/// \brief Read comments on CostFuncDiff::gradient()
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> gradient();


//-----------------------------------------------------------------------
/// \brief The residual vector function 
///
/// This method must be implemented by the classes derived 
/// from this class.
///
/// The parameters (the point in the parameter space) must have
/// already been set before calling this method.  The parameters are
/// already set if one of the following methods is already called
/// successfully:
///   - parameters() (see ProblemState class)
///   - cost_x() (see CostFunc class)
///   - gradient_x() (see CostFuncDiff class)
///   - cost_gradient_x() (see CostFuncDiff class)
///   - residual_x()
///   - jacobian_x()
///   - residual_jacobian_x()
///
/// If the parameters are already set, then this method returns the
/// residual of the NLLS problem at the current set point.
///
/// The size of the residual vector can be obtained in advance
/// by calling residual_size().
/// 
/// \return The residual of the NLLS problem
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> residual() = 0;


//-----------------------------------------------------------------------
/// \brief The residual function with parameters
///
/// This method also evaluates the residual of the NLLS problem;
/// however, it sets the problem at the input new point and then
/// evaluates the residual.
///
/// The size of the residual vector can be obtained in advance
/// by calling residual_size().
/// 
/// \param[in] x
///            New set of parameters
///
/// \return The residual of the cost function
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> residual_x(const blitz::Array<double, 1>& x)
  { parameters(x); return residual(); }


//-----------------------------------------------------------------------
/// \brief The Jacobian matrix function 
///
/// This method must be implemented by the classes derived 
/// from this class.
///
/// The parameters (the point in the parameter space) must have
/// already been set before calling this method.  The parameters are
/// already set if one of the following methods is already called
/// successfully:
///   - parameters() (see ProblemState class)
///   - cost_x() (see CostFunc class)
///   - gradient_x() (see CostFuncDiff class)
///   - cost_gradient_x() (see CostFuncDiff class)
///   - residual_x()
///   - jacobian_x()
///   - residual_jacobian_x()
///
/// If the parameters are already set, then this method returns the
/// Jacobian of the residual of the NLLS problem at the current set
/// point.
///
/// The sizes of the Jacobian matrix can be obtained in advance:
///   - The number of its rows is the same as residual_size().
///   - The number of its columns is the same as 
///     ProblemState::expected_parameter_size() or
///     CostFuncDiff::gradient_size()
/// 
/// \return The Jacobian of the residual of the NLLS problem
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> jacobian() = 0;


//-----------------------------------------------------------------------
/// \brief The Jacobian function with parameters
///
/// This method also evaluates the Jacobian of the residual of the 
/// NLLS problem; however, it sets the problem at the input new point
/// and then evaluates the Jacobian.
///
/// The sizes of the Jacobian matrix can be obtained in advance as
/// mentioned in the comments on jacobian() method.
/// 
/// \param[in] x
///            New set of parameters
///
/// \return The Jacobian of the residual of the NLLS problem
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> jacobian_x(const blitz::Array<double, 1>& x)
  { parameters(x); return jacobian(); }


//-----------------------------------------------------------------------
/// \brief The residual function and its Jacobian together
///
/// This method passes to the caller the evaluated residual function
/// and its Jacobian at the current set point.
///
/// The parameters (the point in the parameter space) must have
/// already been set before calling this method.  The parameters are
/// already set if one of the following methods is already called
/// successfully:
///   - parameters() (see ProblemState class)
///   - cost_x() (see CostFunc class)
///   - gradient_x()
///   - cost_gradient_x()
///   - residual_x()
///   - jacobian_x()
///   - residual_jacobian_x()
///
/// \param[out] r
///             The residual vector
///
/// \param[out] j
///             The Jacobian matrix
//-----------------------------------------------------------------------

  virtual void residual_jacobian(
     blitz::Array<double, 1>& r, blitz::Array<double, 2>& j);


//-----------------------------------------------------------------------
/// \brief The residual and its Jacobian with parameters
///
/// This method passes to the caller the evaluated residual function
/// and its Jacobian after setting the problem at the input new point.
///
/// \param[in] x
///            New set of parameters
/// 
/// \param[out] r
///             The residual vector
///
/// \param[out] j
///             The Jacobian matrix
//-----------------------------------------------------------------------

  virtual void residual_jacobian_x(const blitz::Array<double, 1>& x,
     blitz::Array<double, 1>& r, blitz::Array<double, 2>& j)
  { parameters(x); residual_jacobian(r,j); }


//-----------------------------------------------------------------------
/// \brief Returns the number of the times residual has been evaluated.
///
/// \return The number of the times residual has been evaluated
//-----------------------------------------------------------------------

  virtual int num_residual_evaluations() const
  { return num_cost_evaluations(); }


//-----------------------------------------------------------------------
/// \brief Returns the number of the times Jacobian has been evaluated.
///
/// \return The number of the times Jacobian has been evaluated
//-----------------------------------------------------------------------

  virtual int num_jacobian_evaluations() const
  { return num_der1_evaluations(); }


//-----------------------------------------------------------------------
/// \brief The size of the residual returned by residual()
///
/// This method must be implemented by the classes derived from this
/// class.
///
/// \return The size of the residual that will be returned by residual()
//-----------------------------------------------------------------------

  virtual int residual_size() const = 0;


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSProblem"; }


protected:

};
}
#endif
