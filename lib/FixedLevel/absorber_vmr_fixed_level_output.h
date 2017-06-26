#ifndef ABSORBER_VMR_FIXED_LEVEL_OUTPUT_H
#define ABSORBER_VMR_FIXED_LEVEL_OUTPUT_H
#include "register_output_base.h"
#include "absorber_vmr_fixed_level.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the AbsorberVmrFixedLevel class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the AbsorberVmrFixedLevel class.
*******************************************************************/
class AbsorberVmrFixedLevelOutput : public RegisterOutputBase {
public:
  AbsorberVmrFixedLevelOutput(const boost::shared_ptr<AbsorberVmrFixedLevel>& A,
		      const boost::shared_ptr<StateVector>& Sv)
    : a(A), sv(Sv), num_level(A->pressure()->max_number_level()) {}
  virtual ~AbsorberVmrFixedLevelOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<AbsorberVmrFixedLevel> a;
  boost::shared_ptr<StateVector> sv;
  int num_level;
};
}
#endif
