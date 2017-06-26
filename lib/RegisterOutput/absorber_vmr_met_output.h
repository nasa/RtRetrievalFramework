#ifndef ABSORBER_VMR_MET_OUTPUT_H
#define ABSORBER_VMR_MET_OUTPUT_H
#include "register_output_base.h"
#include "absorber_vmr_met.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the AbsorberVmrMet class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the AbsorberVmrMet class.
*******************************************************************/
class AbsorberVmrMetOutput : public RegisterOutputBase {
public:
  AbsorberVmrMetOutput(const boost::shared_ptr<AbsorberVmrMet>& A)
    : a(A) {}
  virtual ~AbsorberVmrMetOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<AbsorberVmrScaled> a;
};
}
#endif
