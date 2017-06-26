#ifndef ABSORBER_VMR_LEVEL_SCALED_OUTPUT_H
#define ABSORBER_VMR_LEVEL_SCALED_OUTPUT_H
#include "register_output_base.h"
#include "absorber_vmr_level_scaled.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the AbsorberVmrLevelScaled class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the AbsorberVmrLevelScaled class.
*******************************************************************/
class AbsorberVmrLevelScaledOutput : public RegisterOutputBase {
public:
  AbsorberVmrLevelScaledOutput(const boost::shared_ptr<AbsorberVmrLevelScaled>& A)
    : a(A) {}
  virtual ~AbsorberVmrLevelScaledOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<AbsorberVmrScaled> a;
};
}
#endif
