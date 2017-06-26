#ifndef CONNOR_CONVERGENCE_OUTPUT_H
#define CONNOR_CONVERGENCE_OUTPUT_H
#include "register_output_base.h"
#include "connor_convergence.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the ConnorConvergence class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the ConnorConvergence class.
*******************************************************************/
class ConnorConvergenceOutput : public RegisterOutputBase {
public:
  ConnorConvergenceOutput(const boost::shared_ptr<ConnorConvergence>& Conv) 
    : conv(Conv) {}
  virtual ~ConnorConvergenceOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<ConnorConvergence> conv;
};
}
#endif
