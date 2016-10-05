#ifndef PROBLEM_STATE_H
#define PROBLEM_STATE_H
#include <printable.h>
#include <blitz/array.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for all problem states
///
/// ProblemState is the base class for the state of all optimization
/// problems.
///
/// An optimization problem is just a cost function to be minimized. 
/// Therefore, in its simplest form, the state of an optimization 
/// problem is the point in the parameter space where the cost
/// function is currently evaluated.  In other words, given a point
/// in the parameter space, everything else needed (cost function
/// value and its derivatives if needed) can be determined. 
/// Therefore, in the context of our optimization problem we may
/// also use "problem state" to refer to the point in the parameter
/// space where the optimization related cost function is evaluated.
///
/// However, the role of this class is expanded and has resulted in 
/// a class hierarchy rooted at ProblemState class.  It could
/// be very expensive to evaluate some cost functions and/or their
/// derivatives; therefore, it is desirable to store computationally
/// expensive components of a cost function after evaluation.  The
/// classes in the class hierarchy (rooted at ProblemState class) 
/// enable maintaining of the problem state (the current point in the
/// parameter space) and computationally expensive components of the
/// cost function just in case they are needed repeatedly.
///
/// This hierarchy provides a systematic way to 
///   - store an optimization problem state and the related 
///     computationally expensive components,
///   - delete all stored components of the cost function 
///     when the state changes, and
///   - determine when a new state is different enough from
///     the current state to be considered a change in the state.
///
/// All optimization problem classes in the class hierarchy rooted
/// at CostFunc class must directly or indirectly inherit at least
/// ProblemState class.  Given that CostFunc is derived form 
/// ProblemState, then all classes in the problem class hierarchy 
/// will inherit CostFunc automatically.  However, a problem class
/// may also optionally inherit another appropriate problem state
/// class in ProblemState class hierarchy.
//-----------------------------------------------------------------------

class ProblemState : 
    public Printable<ProblemState> {

public:


//-----------------------------------------------------------------------
/// \brief Default constructor
//-----------------------------------------------------------------------

  ProblemState() {}


//-----------------------------------------------------------------------
/// \brief Copy constructor
/// 
/// \param[in] s 
///            another ProblemState
//-----------------------------------------------------------------------

  ProblemState(const ProblemState& s)
  { set(s); }


  virtual ~ProblemState() {}


//-----------------------------------------------------------------------
/// \brief Makes self a copy of the input state
///
/// This method makes the object, for which it is called, a copy of
/// the input state.
/// 
/// \param[in] s 
///            another ProblemState
//-----------------------------------------------------------------------

  virtual void set(const ProblemState& s)
  { X.reference(s.X.copy()); }


//-----------------------------------------------------------------------
/// \brief Deletes data contents.
///
/// This method deletes state.  It must be reimplemented by other
/// classes derived from this class to delete other saved 
/// components associated with the state as well.
//-----------------------------------------------------------------------

  virtual void clear()
  { X.free(); }


//-----------------------------------------------------------------------
/// \brief Checks whether or not new input parameters are different
///        from the current ones
///
/// The methods checks to see whether or not the new input parameters
/// (point in the parameter space) are different from the parameters
/// maintained by the object for which the method is called.
///
/// If the size of the input parameters is not equal to the expected
/// size of the parameters (check comments on expected_parameter_size),
/// then the method will throw an exception.
///
/// If the object for which the method is called has currently no
/// parameters set, then the method returns true.  Otherwise, the 
/// method uses some algorithm to figure out when the difference
/// is "big enough" to be considered different. If the method 
/// determines that the new input parameters are different from the 
/// current parameters, then it returns true, otherwise, it returns
/// false.
///
/// \param[in] x
///            New set of parameters
//-----------------------------------------------------------------------

  virtual bool parameters_different(const blitz::Array<double, 1>& x) const;


//-----------------------------------------------------------------------
/// \brief Sets the problem at a new point in the parameter space.
///
/// The method calls parameters_different() to determine whether or
/// not the new parameters are different:
///   - If different, then it deletes the object state (see clear()),
///     and the input x becomes the current state.
///   - If not different, then the input parameters are ignored.
///
/// \param[in] x 
///            New set of parameters
//-----------------------------------------------------------------------

  virtual void parameters(const blitz::Array<double, 1>& x);


//-----------------------------------------------------------------------
/// \brief Returns the current parameters.
/// 
/// \return Current parameter
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> parameters() const
  { return X.copy(); }


//-----------------------------------------------------------------------
/// \brief Returns the size of the parameters.
///
/// \return Size of parameters
//-----------------------------------------------------------------------

  virtual int parameter_size() const
  { return X.rows(); }


//-----------------------------------------------------------------------
/// \brief Returns the expected size of the parameters.
///
/// This method must be reimplemented by the problem class the
/// inherits ProblemState.  It is only in the context of an
/// optimization problem that one knows what the size of the 
/// parameters (number of the dimensions of the parameter space)
/// is.
///
/// This method is intentionally implemented here instead of being
/// left as a pure virtual method.  The intention is that the user
/// to be able to create an object of this class or its derived 
/// classes for the purpose of preserving an older state of a 
/// problem if needed.
///
/// \return Expected size of parameters
//-----------------------------------------------------------------------

  virtual int expected_parameter_size() const
  { return -1; }


//-----------------------------------------------------------------------
/// \brief Checks that the parameters are set correctly.
///
/// This method checks to see whether or not the parameters are set
/// correctly.  If the parameters are not set correctly then it throws 
/// an exception.
//-----------------------------------------------------------------------

  virtual void assert_parameter_set_correctly() const
  { assert_parameter_correct(X); }


//-----------------------------------------------------------------------
/// \brief Checks that the new input parameters are correct.
///
/// This method checks to see whether or not the new input parameters
/// are correct.  If the parameters are not correct then it throws 
/// an exception.
//-----------------------------------------------------------------------

  virtual void assert_parameter_correct(const blitz::Array<double, 1>& x) const;


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "ProblemState"; }


protected:

  blitz::Array<double, 1> X;

};
}
#endif
