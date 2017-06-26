// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "solver_iteration_log.h"
%}
%base_import(connor_solver)
%import "state_vector.i"
%fp_shared_ptr(FullPhysics::SolverIterationLog);

namespace FullPhysics {
  class SolverIterationLog : public Observer<ConnorSolver> {
  public:
    SolverIterationLog(const boost::shared_ptr<StateVector>& Sv);
    void notify_update(const ConnorSolver& solver);
    std::string print_to_string() const;
};
}
