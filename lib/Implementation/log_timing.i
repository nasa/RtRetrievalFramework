// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "log_timing.h"
%}
%base_import(connor_solver)
%fp_shared_ptr(FullPhysics::LogTiming);

namespace FullPhysics {
class LogTiming : public Observer<ConnorSolver> {
public:
  std::string print_to_string();
  LogTiming();
  virtual void notify_update(const ConnorSolver& Solver);
  virtual void write_to_log(const std::string& Prefix = "") const;
};
}

