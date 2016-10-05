#ifndef ABSORBER_VMR_LEVEL_OUTPUT_H
#define ABSORBER_VMR_LEVEL_OUTPUT_H
#include "register_output_base.h"
#include "absorber_vmr_level.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the AbsorberVmrLevel class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the AbsorberVmrLevel class.
*******************************************************************/
class AbsorberVmrLevelOutput : public RegisterOutputBase {
public:
  AbsorberVmrLevelOutput(const boost::shared_ptr<AbsorberVmrLevel>& A)
    : a(A) {}
  virtual ~AbsorberVmrLevelOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<AbsorberVmrLevel> a;
};
}
#endif
