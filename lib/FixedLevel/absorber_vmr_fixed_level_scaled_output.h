#ifndef ABSORBER_VMR_FIXED_LEVEL_SCALED_OUTPUT_H
#define ABSORBER_VMR_FIXED_LEVEL_SCALED_OUTPUT_H
#include "register_output_base.h"
#include "absorber_vmr_fixed_level_scaled.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the AbsorberVmrFixedLevel class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the AbsorberVmrFixedLevel class.
*******************************************************************/
class AbsorberVmrFixedLevelScaledOutput : public RegisterOutputBase {
public:
  AbsorberVmrFixedLevelScaledOutput(const boost::shared_ptr<AbsorberVmrFixedLevelScaled>& A)
    : a(A), num_level(A->pressure()->max_number_level()) {}
  virtual ~AbsorberVmrFixedLevelScaledOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<AbsorberVmrFixedLevelScaled> a;
  int num_level;
};
}
#endif
