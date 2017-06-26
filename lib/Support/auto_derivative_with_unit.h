#ifndef AUTO_DERIVATIVE_WITH_UNIT_H
#define AUTO_DERIVATIVE_WITH_UNIT_H
#include "printable.h"
#include "auto_derivative.h"
#include "double_with_unit.h"
#include "unit.h"

namespace FullPhysics {
/****************************************************************//**
  This is a AutoDerivative that also has units associated with it.
  This is a simple structure that just keeps these two things
  together. 
*******************************************************************/

template<class T> class AutoDerivativeWithUnit: 
    public Printable<AutoDerivativeWithUnit<T> >,
  boost::field_operators<AutoDerivativeWithUnit<T> >
{
public:
  AutoDerivativeWithUnit() {}
  AutoDerivativeWithUnit(const AutoDerivative<T>& V, const Unit& U)
    : value(V), units(U) {}
  AutoDerivativeWithUnit(const AutoDerivative<T>& V, const std::string& U)
    : value(V), units(U) {}
  AutoDerivativeWithUnit(const AutoDerivative<T>& V)
    : value(V), units(units::dimensionless) {}
  AutoDerivativeWithUnit(const T& V)
    : value(V), units(units::dimensionless) {}
  AutoDerivativeWithUnit(const DoubleWithUnit& V)
    : value(V.value), units(V.units) {}
  AutoDerivative<T> value;
  Unit units;

//-----------------------------------------------------------------------
/// Basic math operators for class.
//-----------------------------------------------------------------------

  inline  AutoDerivativeWithUnit<T>& 
  operator*=(const AutoDerivativeWithUnit<T>& D)
  { value *= D.value; units *= D.units; return *this;}
  inline AutoDerivativeWithUnit<T>& 
  operator/=(const AutoDerivativeWithUnit<T>& D)
  { value /= D.value; units /= D.units; return *this;}
  inline AutoDerivativeWithUnit<T>& 
  operator+=(const AutoDerivativeWithUnit<T>& D)
  { value += D.value * FullPhysics::conversion(D.units, units); return *this;}
  inline AutoDerivativeWithUnit<T>& 
  operator-=(const AutoDerivativeWithUnit<T>& D)
  { value -= D.value * FullPhysics::conversion(D.units, units); return *this;}

//-----------------------------------------------------------------------
/// Convert to the given units.
//-----------------------------------------------------------------------

  inline AutoDerivativeWithUnit<T> convert(const Unit& R) const
  { return AutoDerivativeWithUnit<T>(value * FullPhysics::conversion(units, R), R); }

  void print(std::ostream& Os) const 
  { Os << value << " " << units.name(); }
};
}
#endif
