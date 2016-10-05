#ifndef NLLS_MAX_A_POSTERIORI_H
#define NLLS_MAX_A_POSTERIORI_H
#include <max_a_posteriori.h>
#include <nlls_problem.h>
#include <nlls_problem_state.h>
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/******************************************************************
  This class casts a maximum a posteriori problem into the
  Nonlinear Least Squares problem.
*******************************************************************/
class NLLSMaxAPosteriori: 
    public NLLSProblem, public NLLSProblemState {

public:

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

  NLLSMaxAPosteriori(const boost::shared_ptr<MaxAPosteriori>& map, bool together=false)
    : NLLSProblem(), MAP(map), Compute_together(together)
  {}

  virtual ~NLLSMaxAPosteriori() {}


//-----------------------------------------------------------------------
/// Return the residual of the NLLS problem
/// at the current set point
/// 
/// \return Residual
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> residual();


//-----------------------------------------------------------------------
/// Return the Jacobian of the residual of the NLLS problem
/// at the current set point.
///
/// \return The Jacobian of the cost function.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> jacobian();


//-----------------------------------------------------------------------
/// Return the size of the residual that will be returned by residual()
//-----------------------------------------------------------------------

  virtual int residual_size() const
  { return MAP->measurement_size()+MAP->expected_parameter_size(); }


//-----------------------------------------------------------------------
/// Return the expected size of the parameter X.
//-----------------------------------------------------------------------

  virtual int expected_parameter_size() const
  { return MAP->expected_parameter_size(); }


//-----------------------------------------------------------------------
/// Sets the problem at a new point in the parameter space.
/// 
/// \param x Input value
//-----------------------------------------------------------------------

  virtual void parameters(const blitz::Array<double, 1>& x);


//-----------------------------------------------------------------------
/// Just returns the current values of parameters.
/// This method is redefined here (see the root base
/// class) because of a compiler bug; otherwise, there
/// should be no need for its redefinition.
/// 
/// \return Current parameter values
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> parameters() const
  { return NLLSProblemState::parameters(); }



  boost::shared_ptr<MaxAPosteriori> max_a_posteriori()
  { return MAP; }


//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSMaxAPosteriori"; }


protected:

  const boost::shared_ptr<MaxAPosteriori> MAP;
  bool Compute_together;

};
}
#endif
