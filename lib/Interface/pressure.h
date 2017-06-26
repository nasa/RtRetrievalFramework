#ifndef PRESSURE_H
#define PRESSURE_H
#include "state_vector.h"
#include "observer.h"
#include "array_ad_with_unit.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the pressure portion of the state.

  Note that in a retrieval, there are typically two different pressure
  levels of interest. On is the pressure levels where various initial
  parameters are defined, e.g. Temperature read from an ECMWF file at
  specific pressure levels. The second set is the current pressure
  levels that define the layers used in the Radiative Transfer
  calculation.  The first set is fixed constant level, it is whatever
  was used when we initial read the input data. The second will
  potentially vary as we do a retrieval.

  This class captures the second, potentially varying set of pressure.

  Other classes typically depend on the pressure levels, e.g.,
  Altitude. As a convenience to these classes, the Pressure class can
  notify them when it is changed. These classes can register
  themselves as Observers of the Pressure object if desired.

  When implementing a new class, you almost always will want to derive
  from PressureImpBase rather than from this class. See that class for
  a description.
*******************************************************************/

class Pressure : virtual public StateVectorObserver, 
		 public Observable<Pressure> {
public:
  virtual ~Pressure() {}
  virtual void add_observer(Observer<Pressure>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Pressure>& Obs) 
  { remove_observer_do(Obs, *this);}

  AutoDerivativeWithUnit<double> surface_pressure() const;

//-----------------------------------------------------------------------
/// Return the current surface pressure value, without the
/// gradient. This is in Pascals. 
//-----------------------------------------------------------------------

  double surface_pressure_value() const 
  {return surface_pressure().convert(units::Pa).value.value();}

//-----------------------------------------------------------------------
/// This returns the pressure grid to use for layer retrieval, along
/// with the gradient of each of the pressure grid values with the
/// state vector.
//-----------------------------------------------------------------------

  virtual ArrayAdWithUnit<double, 1> pressure_grid() const = 0;

//-----------------------------------------------------------------------
/// This is the number of layers. This is the same as
/// pressure_grid.rows() - 1.
//-----------------------------------------------------------------------

  int number_layer() const {return pressure_grid().value.rows() - 1;}

//-----------------------------------------------------------------------
/// This is the number of levels. This is the same as
/// pressure_grid.rows() - 1.
//-----------------------------------------------------------------------

  int number_level() const {return pressure_grid().value.rows();}

//-----------------------------------------------------------------------
/// The maximum number of levels that we can have. The default is just
/// number_level() (i.e., we don't change the number of levels from
/// one iteration to the next).
//-----------------------------------------------------------------------

  virtual int max_number_level() const {return number_level();}

//-----------------------------------------------------------------------
/// Clone a Pressure object. Note that the cloned version will *not*
/// be attached to a StateVector or Observer<Pressure>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" Pressure object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Pressure> clone() const = 0;
};
}
#endif
