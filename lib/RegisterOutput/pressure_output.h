#ifndef PRESSURE_OUTPUT_H
#define PRESSURE_OUTPUT_H
#include "register_output_base.h"
#include "pressure.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the Pressure class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the Pressure class.
*******************************************************************/
class PressureOutput : public RegisterOutputBase {
public:
  PressureOutput(const boost::shared_ptr<Pressure>& P,
		 const boost::shared_ptr<StateVector>& Sv)
    : p(P), sv(Sv) {}
  virtual ~PressureOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<Pressure> p;
  boost::shared_ptr<StateVector> sv;
};
}
#endif
