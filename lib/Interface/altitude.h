#ifndef ALTITUDE_H
#define ALTITUDE_H
#include "observer.h"
#include "array_ad.h"
#include "auto_derivative_with_unit.h"

namespace FullPhysics {
class Temperature;
class Pressure;
/****************************************************************//**
   The class handles the calculation of the altitude and gravity
   constants. 

  Other objects may depend on the altitude, and should be updated
  when the altitude is updated. To facilitate that, this class in
  an Oberverable, and objects can add themselves as Observers to be
  notified when the temperature is updated.
*******************************************************************/
class Altitude : public Observable<Altitude>, public Printable<Altitude> {
public:
  virtual ~Altitude() {}
  virtual void add_observer(Observer<Altitude>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Altitude>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Return altitude grid for the given pressure.
//-----------------------------------------------------------------------

  virtual AutoDerivativeWithUnit<double> 
  altitude(const AutoDerivativeWithUnit<double>& P) const = 0;

//-----------------------------------------------------------------------
/// Return gravity constant for the given pressure. 
//-----------------------------------------------------------------------

  virtual AutoDerivativeWithUnit<double> 
  gravity(const AutoDerivativeWithUnit<double>& P) 
    const = 0;

  virtual void print(std::ostream& Os) const { Os << "Altitude";}

//-----------------------------------------------------------------------
/// Clone an Altitude object. Note that the cloned version will *not*
/// be attached to and StateVector or Observer<Altitude>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" Altitude object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Altitude> clone() const = 0;

//-----------------------------------------------------------------------
/// This version of clone takes a pressure and temperature to use. The
/// intent is that the pressure and temperature has been cloned from
/// the original pressure and temperature (although this class has no
/// way to verify this). This allows sets of objects to be cloned
/// using a common Pressure and Temperature clone, e.g. Atmosphere.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Altitude> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<Temperature>& Temp) const = 0;
};
}
#endif
