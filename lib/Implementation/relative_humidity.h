#ifndef RELATIVE_HUMIDITY_H
#define RELATIVE_HUMIDITY_H
#include "absorber.h"
#include "temperature.h"
#include "pressure.h"

namespace FullPhysics {
/****************************************************************//**
   This calculates the relative humidity. This is a straight forward
   calculation from the H2O AbsorberVmr, Temperature, and Pressure. 
   There is no existing class that this obviously goes to, so we 
   just pull this out into its own class.
*******************************************************************/
class RelativeHumidity : public Printable<RelativeHumidity> {
public:
//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

  RelativeHumidity(const boost::shared_ptr<Absorber>& Abs, 
		   const boost::shared_ptr<Temperature>& Temp,
		   const boost::shared_ptr<Pressure>& Press)
    : absorber(Abs), temp(Temp), press(Press)
  { }
  virtual ~RelativeHumidity() {}
  void print(std::ostream& Os);
  ArrayAd<double, 1> relative_humidity_grid() const;
private:
  boost::shared_ptr<Absorber> absorber;
  boost::shared_ptr<Temperature> temp;
  boost::shared_ptr<Pressure> press;
};
}
#endif
