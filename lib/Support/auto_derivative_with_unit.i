// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "auto_derivative_with_unit.h"
%}

%base_import(generic_object)
%import "unit.i"
%import "auto_derivative.i"

namespace FullPhysics {
template<class T> class AutoDerivativeWithUnit : public GenericObject
{
public:
  AutoDerivativeWithUnit();
  AutoDerivativeWithUnit(const AutoDerivative<T>& V, const Unit& U);
  AutoDerivativeWithUnit(const AutoDerivative<T>& V, const std::string& U);
  AutoDerivativeWithUnit(const AutoDerivative<T>& V);
  AutoDerivativeWithUnit<T> convert(const Unit& R) const;
  AutoDerivativeWithUnit<T> convert(const std::string& R) const;
  std::string print_to_string() const;
  %extend {
    FullPhysics::AutoDerivative<T> _value() const { return $self->value;}
    void _value_set(const FullPhysics::AutoDerivative<T>& V) { $self->value = V;}
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
}
%fp_shared_ptr(FullPhysics::AutoDerivativeWithUnit<double>)
%template(AutoDerivativeWithUnitDouble) FullPhysics::AutoDerivativeWithUnit<double>;




