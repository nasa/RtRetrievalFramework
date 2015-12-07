#ifndef TEMPERATURE_H
#define TEMPERATURE_H
#include "state_vector.h"
#include "observer.h"
#include "auto_derivative.h"
#include "array_with_unit.h"
#include "pressure.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the temperature portion of the state.

  Other objects may depend on the temperature, and should be updated
  when the temperature is updated. To facilitate that, this class in
  an Oberverable, and objects can add themselves as Observers to be
  notified when the temperature is updated.

  When implementing a new class, you almost always will want to derive
  from TemperatureImpBase rather than from this class. See that class for
  a description.
*******************************************************************/
class Temperature : virtual public StateVectorObserver,
		    public Observable<Temperature> {
public:
  virtual ~Temperature() {}
  virtual void add_observer(Observer<Temperature>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Temperature>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// The temperature can vary quickly over a small pressure range,
/// e.g. at the tropopause and stratopause. It is important that this
/// structure is included in anything using the temperature, e.g., the
/// integration does to calculate the optical depth of a layer in
/// AbsorberAbsco. 
///
/// This supplied "important" pressures where something interesting in
/// the temperature may be happening. 
///
/// The default is that there are not important pressures, but a
/// derived class can override this, e.g. give the ECMWF pressure
/// levels.
//-----------------------------------------------------------------------

  virtual ArrayWithUnit<double, 1> important_pressure_level() const 
  {
    blitz::Array<double, 1> empty;
    return ArrayWithUnit<double, 1>(empty, units::Pa);
  }

//-----------------------------------------------------------------------
/// Return the temperature at the given pressure (in Pascals)
///
/// This is in Kelvin.
//-----------------------------------------------------------------------

  virtual AutoDerivativeWithUnit<double> 
  temperature(const AutoDerivativeWithUnit<double>& Press) const = 0;

  virtual ArrayAdWithUnit<double, 1> temperature_grid(const Pressure& P) const;

//-----------------------------------------------------------------------
/// Clone a Temperature object. Note that the cloned version will *not*
/// be attached to a StateVector or Observer<Temperature>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" Temperature object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Temperature> clone() const = 0;

//-----------------------------------------------------------------------
/// This version of clone takes a pressure to use. The intent is that
/// the pressure has been cloned from the original pressure (although
/// this class has no way to verify this). This allows sets of objects
/// to be cloned using a common Pressure clone, e.g. Atmosphere.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Temperature> 
  clone(const boost::shared_ptr<Pressure>& Press) const = 0;
};
}
#endif
