#ifndef TEMPERATURE_FIXED_LEVEL_OUTPUT_H
#define TEMPERATURE_FIXED_LEVEL_OUTPUT_H
#include "register_output_base.h"
#include "temperature_fixed_level.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the TemperatureFixedLevel class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the TemperatureFixedLevel class.
*******************************************************************/
class TemperatureFixedLevelOutput : public RegisterOutputBase {
public:
  TemperatureFixedLevelOutput(const boost::shared_ptr<TemperatureFixedLevel>& T) 
  : t(T) {}
  virtual ~TemperatureFixedLevelOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<TemperatureFixedLevel> t;
};
}
#endif
