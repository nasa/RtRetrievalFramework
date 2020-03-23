#ifndef ABSORBER_VMR_LOG_LEVEL_OUTPUT_H
#define ABSORBER_VMR_LOG_LEVEL_OUTPUT_H
#include "register_output_base.h"
#include "absorber_vmr_log_level.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the AbsorberVmrLogLevel class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the AbsorberVmrLogLevel class.
*******************************************************************/
class AbsorberVmrLogLevelOutput : public RegisterOutputBase {
public:
  AbsorberVmrLogLevelOutput(const boost::shared_ptr<AbsorberVmrLogLevel>& A)
    : a(A) {}
  virtual ~AbsorberVmrLogLevelOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<AbsorberVmrLogLevel> a;
};
}
#endif
