// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "double_with_unit.h"
%}

%base_import(generic_object)
%import "unit.i"
%fp_shared_ptr(FullPhysics::DoubleWithUnit)
namespace FullPhysics {
class DoubleWithUnit : public GenericObject {
public:
  std::string print_to_string() const;
  DoubleWithUnit();
  DoubleWithUnit(double V, const Unit& U);
  DoubleWithUnit(double V, const std::string& U);
  DoubleWithUnit(double V);
  DoubleWithUnit& operator*=(const DoubleWithUnit& D);
  DoubleWithUnit& operator/=(const DoubleWithUnit& D);
  DoubleWithUnit& operator+=(const DoubleWithUnit& D);
  DoubleWithUnit& operator-=(const DoubleWithUnit& D);
  DoubleWithUnit convert(const Unit& R) const;
  DoubleWithUnit convert(const std::string& R) const;
  DoubleWithUnit convert_wave(const Unit& R) const;
  DoubleWithUnit convert_wave(const std::string& R) const;
  // There is a cyclic dependency here. Go ahead and just leave this
  // out for now, we can come up with a solution if we need to.
  // DoubleWithUnit convert_wave(const Unit& R, const SpectralDomain& Sd) const;
  %extend {
    double _value() const {return $self->value; }
    void _value_set(double V) { $self->value = V;}
    Unit _units() const {return $self->units; }
    void _units_set(const FullPhysics::Unit& U) {$self->units = U;}
    DoubleWithUnit __mul__(const DoubleWithUnit& Y) 
    { return *$self * Y; }
    // Python 2 division operator name
    DoubleWithUnit __div__(const DoubleWithUnit& Y) 
    { return *$self / Y; }
    // Python 3 division operator name
    DoubleWithUnit __truediv__(const DoubleWithUnit& Y) 
    { return *$self / Y; }
    DoubleWithUnit __add__(const DoubleWithUnit& Y) 
    { return *$self + Y; }
    DoubleWithUnit __sub__(const DoubleWithUnit& Y) 
    { return *$self - Y; }
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

