// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "connor_solver_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "connor_solver.i"

%fp_shared_ptr(FullPhysics::ConnorSolverOutput);

namespace FullPhysics {
class ConnorSolverOutput : public RegisterOutputBase {
public:
  ConnorSolverOutput(const boost::shared_ptr<ConnorSolver>& Solver, 
		     bool Write_jacobian = false);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
};
}


