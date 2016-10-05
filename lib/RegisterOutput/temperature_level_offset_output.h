#ifndef TEMPERATURE_LEVEL_OFFSET_OUTPUT_H
#define TEMPERATURE_LEVEL_OFFSET_OUTPUT_H
#include "register_output_base.h"
#include "temperature_level_offset.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the TemperatureLevelOffset class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the TemperatureLevelOffset class.
*******************************************************************/
class TemperatureLevelOffsetOutput : public RegisterOutputBase {
public:
  TemperatureLevelOffsetOutput(const boost::shared_ptr<TemperatureLevelOffset>& T) 
  : t(T) {}
  virtual ~TemperatureLevelOffsetOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<TemperatureOffset> t;
};
}
#endif
