// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "connor_convergence_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "connor_convergence.i"

%fp_shared_ptr(FullPhysics::ConnorConvergenceOutput);

namespace FullPhysics {
class ConnorConvergenceOutput : public RegisterOutputBase {
public:
  ConnorConvergenceOutput(const boost::shared_ptr<ConnorConvergence>& Conv);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
};
}


