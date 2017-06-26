#ifndef COST_FUNC_DIFF_H
#define COST_FUNC_DIFF_H
#include <cost_func.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for all problem classes that implement
///        a cost function and its gradient.
///
/// The class CostFuncDiff is the base class for all problem classes
/// that implement a cost function (a scalar real function with a
/// range that never includes negative numbers) and its gradient
/// (first order derivatives).
///
/// A cost function for which we can compute its gradient can be
/// optimized by methods such as the steepest-descent-method and
/// the conjugate-gradient-method as well as other optimization
/// methods that only use the cost function (no derivatives).
//-----------------------------------------------------------------------

class CostFuncDiff : public CostFunc {

public:


//-----------------------------------------------------------------------
/// \brief Default Constructor
//-----------------------------------------------------------------------

  CostFuncDiff() 
    : CostFunc(), d_count(0)
  {}


  virtual ~CostFuncDiff() {}


//-----------------------------------------------------------------------
/// \brief The gradient of the cost function
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
///   - gradient_x()
///   - cost_gradient_x()
///
/// If the parameters are already set, then this method returns the
/// gradient of the cost function at the current set point.
///
/// The size of the gradient vector can be obtained in advance
/// by calling gradient_size().
///
/// \return The gradient of the cost function
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> gradient() = 0;


//-----------------------------------------------------------------------
/// \brief The gradient function with parameters
///
/// This method also evaluates the gradient of the cost function;
/// however, it sets the problem at the input new point and then
/// evaluates the gradient.
///
/// The size of the gradient vector can be obtained in advance
/// by calling gradient_size().
/// 
/// \param[in] x
///            New set of parameters
///
/// \return The gradient of the cost function
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> gradient_x(const blitz::Array<double, 1>& x)
  { parameters(x); return gradient(); }


//-----------------------------------------------------------------------
/// \brief The cost function and its gradient together
///
/// This method passes to the caller the evaluated cost function and
/// its gradient at the current set point.
///
/// The parameters (the point in the parameter space) must have
/// already been set before calling this method.  The parameters are
/// already set if one of the following methods is already called
/// successfully:
///   - parameters() (see ProblemState class)
///   - cost_x() (see CostFunc class)
///   - gradient_x()
///   - cost_gradient_x()
///
/// \param[out] c
///             The cost function value
///
/// \param[out] g
///             The gradient vector with size gradient_size()
//-----------------------------------------------------------------------

  virtual void cost_gradient(
     double& c, blitz::Array<double, 1>& g);


//-----------------------------------------------------------------------
/// \brief The cost function and its gradient with parameters
///
/// This method passes to the caller the evaluated cost function and
/// its gradient after setting the problem at the input new point.
/// 
/// \param[in] x
///            New set of parameters
///
/// \param[out] c
///             The cost function value
///
/// \param[out] g
///             The gradient vector with size gradient_size()
//-----------------------------------------------------------------------

  virtual void cost_gradient_x(const blitz::Array<double, 1>& x,
     double& c, blitz::Array<double, 1>& g)
  { parameters(x); cost_gradient(c,g); }


//-----------------------------------------------------------------------
/// \brief Returns the number of the times gradient has been evaluated.
///
/// \return The number of the times gradient has been evaluated.
//-----------------------------------------------------------------------

  virtual int num_der1_evaluations() const
  { return d_count; }


//-----------------------------------------------------------------------
/// \brief Sets cost and gradient evaluation counters to zero.
//-----------------------------------------------------------------------

  virtual void zero_num_evaluations()
  { CostFunc::zero_num_evaluations(); d_count = 0; }


//-----------------------------------------------------------------------
/// \brief Returns the size of the gradient vector.
///
/// \return The size of the gradient vecotor
//-----------------------------------------------------------------------

  virtual int gradient_size() const
  { return expected_parameter_size(); }


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "CostFuncDiff"; }


protected:


//-----------------------------------------------------------------------
/// \brief Increments (by 1) gradient evaluation counter.
///
/// The developer of a derived class, where the true gradient
/// evaluation is implemented, must also call this method when the
/// gradient (first order derivatives) is evaluated.
///
/// It is possible to call this method in gradient_x() method; however,
/// at this level we do not know when the gradient gets truly
/// evaluated.
///
/// The developer of a derived class can save the expensive components
/// of the gradient evaluation and reuse them when gradient is
/// inquired and the parameters have not changed.  Therefore, 
/// she has the options of calling increment_num_der1_evaluations()
///   -# any time gradient is inquired or
///   -# only when gradient is truly evaluated.
///
/// If X1 and X2 are two different sets of parameters, then after the
/// sequence of method calls
///   - gradient_x(X1)
///   - gradient_x(X1)
///   - gradient_x(X2)
///   - gradient_x(X2)
///
/// with the first design option mentioned above
/// num_der1_evaluations() will return 4, but with the second option
/// mentioned above num_der1_evaluations() will return 2.
///
/// Therefore, in my opinion it was best not to decide when to call
/// this method at this level in CostFuncDiff class.
//-----------------------------------------------------------------------

  virtual void increment_num_der1_evaluations()
  { d_count++; }


//-----------------------------------------------------------------------
/// \brief Sets the gradient evaluation counter to a desired value.
///
/// This method sets the gradient (first order derivatives) evaluation
/// counter to a desired value.  It is just provided if the developers
/// of the derived classes have any use for it
///
/// \param[in] count
///            Desired value for gradient evaluation counter setting
//-----------------------------------------------------------------------

  virtual void set_num_der1_evaluations(int count)
  { d_count = count; }


private:

  int d_count;

};
}
#endif
