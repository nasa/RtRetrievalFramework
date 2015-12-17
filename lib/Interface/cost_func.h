#ifndef COST_FUNC_H
#define COST_FUNC_H
#include <problem_state.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for all problem classes that implement
///        a cost function.
///
/// The class CostFunc is the base class for all problem classes
/// that implement a cost function (a scalar real function with a
/// range that never includes negative numbers).
///
/// A cost function only (without derivatives of any order) can
/// be optimized by methods such as the Sampling Method that do 
/// not require derivatives of any order.
//-----------------------------------------------------------------------

class CostFunc : 
    virtual public ProblemState {

public:


//-----------------------------------------------------------------------
/// \brief Default Constructor
//-----------------------------------------------------------------------

  CostFunc()
    : c_count(0)
  {}


  virtual ~CostFunc() {}


//-----------------------------------------------------------------------
/// \brief The cost function 
///
/// This method must be implemented by the classes derived 
/// from this class.
///
/// The parameters (the point in the parameter space) must have
/// already been set before calling this method.  The parameters are
/// already set if parameters() method (see ProblemState class) or
/// cost_x() method are already called successfully.  If the 
/// parameters are already set, then this method returns the value of
/// the cost function at the current set point.
/// 
/// \return Cost function value
//-----------------------------------------------------------------------

  virtual double cost() = 0;


//-----------------------------------------------------------------------
/// \brief The cost function with parameters
///
/// This method also evaluates the cost function; however, it sets
/// the problem at the input new point and then evaluates the cost
/// function.
/// 
/// \param[in] x
///            New set of parameters
///
/// \return Cost function value
//-----------------------------------------------------------------------

  virtual double cost_x(const blitz::Array<double, 1>& x)
  { parameters(x); return cost(); }


//-----------------------------------------------------------------------
/// \brief Returns the number of the times cost has been evaluated.
///
/// \return The number of the times cost has been evaluated.
//-----------------------------------------------------------------------

  virtual int num_cost_evaluations() const
  { return c_count; }


//-----------------------------------------------------------------------
/// \brief Sets cost evaluation counter to zero.
//-----------------------------------------------------------------------

  virtual void zero_num_evaluations()
  { c_count = 0; }


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "CostFunc"; }


protected:


//-----------------------------------------------------------------------
/// \brief Increments (by 1) cost evaluation counter.
///
/// The developer of a derived class, where the true cost function
/// evaluation is implemented, must also call this method when the
/// cost function is evaluated.
///
/// It is possible to call this method in cost_x() method; however,
/// at this level we do not know when the cost function gets truly
/// evaluated.
///
/// The developer of a derived class can save the expensive
/// components of a cost function evaluation and reuse them when
/// cost is inquired and the parameters have not changed.  Therefore, 
/// she has the options of calling increment_num_cost_evaluations()
///   -# any time cost is inquired or
///   -# only when cost is truly evaluated.
///
/// If X1 and X2 are two different sets of parameters, then after the
/// sequence of method calls
///   - cost_x(X1)
///   - cost_x(X1)
///   - cost_x(X2)
///   - cost_x(X2)
///
/// with the first design option mentioned above
/// num_cost_evaluations() will return 4, but with the second option
/// mentioned above num_cost_evaluations() will return 2.
///
/// Therefore, in my opinion it was best not to decide when to call
/// this method at this level in CostFunc class.
//-----------------------------------------------------------------------

  virtual void increment_num_cost_evaluations()
  { c_count++; }


//-----------------------------------------------------------------------
/// \brief Sets the cost evaluation counter to a desired value.
///
/// This method sets the cost evaluation counter to a desired value.
/// It is just provided if the developers of the derived classes have
/// any use for it
///
/// \param[in] count
///            Desired value for cost evaluation counter setting
//-----------------------------------------------------------------------

  virtual void set_num_cost_evaluations(int count)
  { c_count = count; }


private:

  int c_count;

};
}
#endif
