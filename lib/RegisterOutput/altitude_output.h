#ifndef ALTITUDE_OUTPUT_H
#define ALTITUDE_OUTPUT_H
#include "register_output_base.h"
#include "altitude.h"
#include "pressure.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the Altitude class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the Pressure class.
*******************************************************************/
class AltitudeOutput : public RegisterOutputBase {
public:
  AltitudeOutput(const boost::shared_ptr<Altitude>& Alt,
                 const boost::shared_ptr<Pressure>& Pres)
    : alt(Alt), pres(Pres) {}
  virtual ~AltitudeOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<Altitude> alt;
  boost::shared_ptr<Pressure> pres;
};
}
#endif
