#ifndef COST_FUNC_STATE_H
#define COST_FUNC_STATE_H
#include <problem_state.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The state for a problem that only its cost function 
///        is implemented
///
/// CostFuncState is used for the problems that only their cost 
/// functions are implemented.  With this class one can store the
/// current point in the parameter space (the state) and the 
/// value of the cost function at that point.
//-----------------------------------------------------------------------

class CostFuncState: 
    virtual public ProblemState {

public:


//-----------------------------------------------------------------------
/// \brief Default constructor
//-----------------------------------------------------------------------

  CostFuncState() {}


//-----------------------------------------------------------------------
/// \brief Copy constructor
/// 
/// \param[in] s 
///            another CostFuncState
//-----------------------------------------------------------------------

  CostFuncState(const CostFuncState& s)
  { set(s); }


  virtual ~CostFuncState() {}


//-----------------------------------------------------------------------
/// \brief Makes self a copy of the input state
///
/// This method makes the object, for which it is called, a copy of
/// the input state.
/// 
/// \param[in] s 
///            another CostFuncState
//-----------------------------------------------------------------------

  virtual void set(const CostFuncState& s);


//-----------------------------------------------------------------------
/// \brief Deletes data contents.
///
/// This method deletes state.  If needed, it must be reimplemented
/// by other classes derived from this class to delete other saved 
/// components associated with the state as well.
//-----------------------------------------------------------------------

  virtual void clear()
  { ProblemState::clear(); C.free(); }


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "CostFuncState"; }


protected:

//-----------------------------------------------------------------------
/// A cost function is a scalar function; therefore, its value at
/// a point is also a scalar.  However, here we use a 1D  array to 
/// store the value of the cost function evaluated at the current 
/// point in the parameter space.  The reason for using a 1D array is
/// only for convenience.  However, the details of how we store the
/// cost function value does not appear in the class interface.
///
/// This must be used as a one element array or an empty array.
//-----------------------------------------------------------------------
  blitz::Array<double, 1> C;

};
}
#endif
