#ifndef NLLS_PROBLEM_SCALEDH
#define NLLS_PROBLEM_SCALEDH
#include <nlls_problem.h>
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/******************************************************************
  The main function of this class is to scale the parameter
  space of a NLLS problem to be solved.  Currently only a
  diagonal scaling matrix is supported.
*******************************************************************/
class NLLSProblemScaled : 
    public NLLSProblem {

public:

//-----------------------------------------------------------------------
/// Default Constructor
//-----------------------------------------------------------------------

  NLLSProblemScaled( const blitz::Array<double, 1>& s,
                     const boost::shared_ptr<NLLSProblem>& p );

  virtual ~NLLSProblemScaled() {}

//-----------------------------------------------------------------------
/// Return the residual of the scaled NLLS problem
/// at the current set point
/// 
/// \return Residual
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> residual();

//-----------------------------------------------------------------------
/// Return the Jacobian of the residual of the scaled
/// NLLS problem at the current set point.
///
/// \return The Jacobian of the residual function.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> jacobian();


  virtual int residual_size() const
  { return P->residual_size(); }


//-----------------------------------------------------------------------
/// Return the size of the parameter X.
//-----------------------------------------------------------------------

  virtual int expected_parameter_size() const
  { return P->expected_parameter_size(); }


//-----------------------------------------------------------------------
/// If x is the input to the NLLS problem that this class 
/// is trying to scale, then this method scales the input
/// to be used by this class, i.e. this->parameters(x).
/// The reason for scaling x outside of this->parameters(x)
/// is that we can also scale an already scaled NLLS problem.
///
/// In summary, the input x is a correct input directly to 
/// the NLLS problem being scaled.  The returned value is
/// correctly scaled to be used as input to this scaled 
/// NLLS problem.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> scale_parameters(const blitz::Array<double, 1>& x) const;


//-----------------------------------------------------------------------
/// The input is correctly scaled to be used as input to
/// this scaled  NLLS problem.  The returned value is a
/// correct direct input to the NLLS problem being scaled.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> unscale_parameters(const blitz::Array<double, 1>& x) const;


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
  { return NLLSProblem::parameters(); }


  boost::shared_ptr<NLLSProblem> nlls_problem()
  { return P; }


//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSProblemScaled"; }

protected:

  blitz::Array<double, 1> S;
  boost::shared_ptr<NLLSProblem> P;

};
}
#endif
