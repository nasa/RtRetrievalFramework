#ifndef PRESSURE_FIXED_LEVEL_OUTPUT_H
#define PRESSURE_FIXED_LEVEL_OUTPUT_H
#include "register_output_base.h"
#include "pressure_fixed_level.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the PressureFixedLevel class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the PressureFixedLevel class.
*******************************************************************/
class PressureFixedLevelOutput : public RegisterOutputBase {
public:
  PressureFixedLevelOutput(const boost::shared_ptr<PressureFixedLevel>& P,
			   const boost::shared_ptr<StateVector>& Sv)
    : p(P), sv(Sv) {}
  virtual ~PressureFixedLevelOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<PressureFixedLevel> p;
  boost::shared_ptr<StateVector> sv;
};
}
#endif
