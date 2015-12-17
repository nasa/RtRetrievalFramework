#ifndef COST_FUNC_DIFF_STATE_H
#define COST_FUNC_DIFF_STATE_H
#include <cost_func_state.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The state for a problem with implemented cost function 
///        and its gradient
///
/// CostFuncDiffState is used for the problems that their cost 
/// functions and their gradients are implemented.  With this class
/// one can store the current point in the parameter space (the state),
/// the value of the cost function at that point, and the gradient of
/// the cost function at the same point.
//-----------------------------------------------------------------------

class CostFuncDiffState : 
    public CostFuncState {

public:


//-----------------------------------------------------------------------
/// \brief Default constructor
//-----------------------------------------------------------------------

  CostFuncDiffState() {}


//-----------------------------------------------------------------------
/// \brief Copy constructor
/// 
/// \param[in] s 
///            another CostFuncDiffState
//-----------------------------------------------------------------------

  CostFuncDiffState(const CostFuncDiffState& s)
  { set(s); }


  virtual ~CostFuncDiffState() {}


//-----------------------------------------------------------------------
/// \brief Makes self a copy of the input state
///
/// This method makes the object, for which it is called, a copy of
/// the input state.
/// 
/// \param[in] s 
///            another CostFuncDiffState
//-----------------------------------------------------------------------

  virtual void set(const CostFuncDiffState& s);


//-----------------------------------------------------------------------
/// \brief Deletes data contents.
///
/// This method deletes state.  If needed, it must be reimplemented
/// by other classes derived from this class to delete other saved 
/// components associated with the state as well.
//-----------------------------------------------------------------------

  virtual void clear()
  { CostFuncState::clear(); G.free(); }


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "CostFuncDiffState"; }


protected:

  blitz::Array<double, 1> G;

};
}
#endif
