#ifndef TEMPERATURE_ECMWF_OUTPUT_H
#define TEMPERATURE_ECMWF_OUTPUT_H
#include "register_output_base.h"
#include "temperature_ecmwf.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the TemperatureEcmwf class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the TemperatureEcmwf class.
*******************************************************************/
class TemperatureEcmwfOutput : public RegisterOutputBase {
public:
  TemperatureEcmwfOutput(const boost::shared_ptr<TemperatureEcmwf>& T) 
  : t(T) {}
  virtual ~TemperatureEcmwfOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<TemperatureOffset> t;
};
}
#endif
