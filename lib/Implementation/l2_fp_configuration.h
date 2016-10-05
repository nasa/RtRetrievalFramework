#ifndef L2_FP_CONFIGURATION_H
#define L2_FP_CONFIGURATION_H
#include "printable.h"
#include "forward_model.h"
#include "initial_guess.h"
#include "logger.h"
#include "initial_guess.h"
#include "connor_solver.h"
#include "output.h"
#include "iterative_solver.h"
#include "max_a_posteriori.h"

namespace FullPhysics {
/****************************************************************//**
   Before running L2 full physics, we need to create the solver that
   we will be using, along with registering whatever output we will
   be generating.

   This class gives the minimum interface needed for the
   configuration, so we can use different methods of actually doing
   this. 
*******************************************************************/

class L2FpConfiguration : public Printable<L2FpConfiguration> {
public:
  virtual ~L2FpConfiguration() {}

//-----------------------------------------------------------------------
/// Logger to use.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<LogImp> logger() const = 0;

//-----------------------------------------------------------------------
/// Forward model. Everything should be initialized to the initial guess.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<ForwardModel> forward_model() const = 0;

//-----------------------------------------------------------------------
/// Solver.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<ConnorSolver> solver() const = 0;

//-----------------------------------------------------------------------
/// Iterative solver.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<IterativeSolver> iterative_solver() const = 0;

//-----------------------------------------------------------------------
/// Maximum a posteriori.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<MaxAPosteriori> max_a_posteriori() const = 0;

//-----------------------------------------------------------------------
/// Initial guess.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<InitialGuess> initial_guess() const = 0;

//-----------------------------------------------------------------------
/// Create output, for both a normal run and for an error run (either
/// or both can be null if we don't want output). This should have all
/// the RegisterOutputBase applied to it that the configuration says
/// should be.
//-----------------------------------------------------------------------

  virtual void output(boost::shared_ptr<Output>& Regular_output,
		      boost::shared_ptr<Output>& Error_output) const = 0;

  virtual void print(std::ostream& Os) const {Os << "L2FpConfiguration";}
};
}
#endif
