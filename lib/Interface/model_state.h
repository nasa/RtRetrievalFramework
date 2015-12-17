#ifndef MODEL_STATE_H
#define MODEL_STATE_H
#include <problem_state.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The state for a parametrized mathematical model
///        (a vector function) and its Jacobian
///
/// This class is used as a state for a mathematical model and not an
/// optimization problem.
///
/// Given NLLSProblemState class it appears thath ModelState is redundant.
/// After all, both classes just are designed to maintain a vector
/// function and its Jacobian.  Then why do we have two similar
/// classes that are only different in their names or the names of
/// some of their members?
///
/// A parametrized mathematical model is not an optimization problem
/// by itself; however, it is a component of an optimization problem
/// when we try to fit the model to measured data.
///
/// When a parametrized mathematical model appears in an optimization
/// problem in the form of a Nonlinear Least Squares problem, the
/// vector model function and its Jacobian are not the same as the
/// vector residual function of the NLLS problem and its Jacobian.
/// They are different, and they are very different when we use 
/// some statistical analysis method to fit the model to the measured
/// data.  In my judgment, emphasizing the 
/// differences and avoiding confusion are more important than 
/// redundancy in this case; therefore, I implemented ModelState
/// as well as NLLSProblemState.
//-----------------------------------------------------------------------

class ModelState : 
    virtual public ProblemState {

public:


//-----------------------------------------------------------------------
/// \brief Default constructor
//-----------------------------------------------------------------------

  ModelState() {}


//-----------------------------------------------------------------------
/// \brief Copy constructor
/// 
/// \param[in] s 
///            another ModelState
//-----------------------------------------------------------------------

  ModelState(const ModelState& s)
  { set(s); }


  virtual ~ModelState() {}


//-----------------------------------------------------------------------
/// \brief Makes self a copy of the input state
///
/// This method makes the object, for which it is called, a copy of
/// the input state.
/// 
/// \param[in] s 
///            another ModelState
//-----------------------------------------------------------------------

  virtual void set(const ModelState& s);


//-----------------------------------------------------------------------
/// \brief Deletes data contents.
///
/// This method deletes state.  If needed, it must be reimplemented
/// by other classes derived from this class to delete other saved 
/// components associated with the state as well.
//-----------------------------------------------------------------------

  virtual void clear()
  { ProblemState::clear(); M.free(); K.free(); }


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "ModelState"; }


protected:

  blitz::Array<double, 1> M;
  blitz::Array<double, 2> K;

};
}
#endif
