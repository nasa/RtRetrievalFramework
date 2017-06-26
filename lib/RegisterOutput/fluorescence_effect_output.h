#ifndef FLUORESCENCE_EFFECT_OUTPUT_H
#define FLUORESCENCE_EFFECT_OUTPUT_H
#include "register_output_base.h"
#include "fluorescence_effect.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the FluorescenceEffect class that
  should be written as output.
*******************************************************************/
class FluorescenceEffectOutput : public RegisterOutputBase {
public:
  FluorescenceEffectOutput(const boost::shared_ptr<FluorescenceEffect>& Fluor)
    : fluor(Fluor) {}
  virtual ~FluorescenceEffectOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<FluorescenceEffect> fluor;
};
}
#endif
