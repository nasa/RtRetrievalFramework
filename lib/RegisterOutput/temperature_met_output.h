#ifndef TEMPERATURE_MET_OUTPUT_H
#define TEMPERATURE_MET_OUTPUT_H
#include "register_output_base.h"
#include "temperature_met.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the TemperatureMet class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the TemperatureMet class.
*******************************************************************/
class TemperatureMetOutput : public RegisterOutputBase {
public:
  TemperatureMetOutput(const boost::shared_ptr<TemperatureMet>& T) 
  : t(T) {}
  virtual ~TemperatureMetOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<TemperatureOffset> t;
};
}
#endif
