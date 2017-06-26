// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "array_ad_with_unit.h"
%}

%import "array_ad.i"
%import "unit.i"

namespace FullPhysics {
template<class T, int D> class ArrayAdWithUnit
{
public:
  ArrayAdWithUnit();
  ArrayAdWithUnit(const ArrayAd<T, D>& V, const Unit& U);
  ArrayAdWithUnit(const ArrayAd<T, D>& V, const std::string& U);
  ArrayAdWithUnit(const ArrayAd<T, D>& V);
  ArrayAdWithUnit<T, D> convert(const Unit& R) const;
  ArrayAdWithUnit<T, D> convert(const std::string& R) const;
  %python_attribute(rows, int)
  %python_attribute(cols, int)
  %python_attribute(depth, int)
  %python_attribute(is_constant, bool)
  %python_attribute(number_variable, int)
  void reference(const ArrayAdWithUnit<T, D>& V);
  %extend {
    FullPhysics::ArrayAd<T, D> _value() const { return $self->value;}
    void _value_set(const FullPhysics::ArrayAd<T, D>& V) { $self->value = V;}
    FullPhysics::Unit _units() const {return $self->units;}
    void _units_set(const FullPhysics::Unit& U) {$self->units = U;}
  }
  %pythoncode {
@property
def value(self):
  return self._value()

@value.setter
def value(self, val):
  self._value_set(val)

@property
def units(self):
  return self._units()

@units.setter
def units(self,val):
    self._units_set(val)
  }
};
%template(ArrayAdWithUnitDouble_1) FullPhysics::ArrayAdWithUnit<double, 1>;
%template(ArrayAdWithUnitDouble_2) FullPhysics::ArrayAdWithUnit<double, 2>;
%template(ArrayAdWithUnitDouble_3) FullPhysics::ArrayAdWithUnit<double, 3>;
%template(ArrayAdWithUnitDouble_4) FullPhysics::ArrayAdWithUnit<double, 4>;
}




