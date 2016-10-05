#ifndef CONNOR_SOLVER_OUTPUT_H
#define CONNOR_SOLVER_OUTPUT_H
#include "register_output_base.h"
#include "connor_solver.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the ConnorSolver class that should be
  written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the ConnnorSolver class.
*******************************************************************/
class ConnorSolverOutput : public RegisterOutputBase {
public:
//-----------------------------------------------------------------------
/// You can optionally turn on the writing of the Jacobian. We don't
/// do this by default because the files get large if we do, and most
/// of the time we don't need the Jacobian values.
//-----------------------------------------------------------------------

  ConnorSolverOutput(const boost::shared_ptr<ConnorSolver>& Solver, 
		     bool Write_jacobian = false)
    : solver(Solver), write_jacobian(Write_jacobian) {}
  virtual ~ConnorSolverOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<ConnorSolver> solver;
  bool write_jacobian;
};
}
#endif
