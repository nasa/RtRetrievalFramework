#ifndef ABSORBER_VMR_ECMWF_OUTPUT_H
#define ABSORBER_VMR_ECMWF_OUTPUT_H
#include "register_output_base.h"
#include "absorber_vmr_ecmwf.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the AbsorberVmrEcmwf class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the AbsorberVmrEcmwf class.
*******************************************************************/
class AbsorberVmrEcmwfOutput : public RegisterOutputBase {
public:
  AbsorberVmrEcmwfOutput(const boost::shared_ptr<AbsorberVmrEcmwf>& A)
    : a(A) {}
  virtual ~AbsorberVmrEcmwfOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<AbsorberVmrScaled> a;
};
}
#endif
