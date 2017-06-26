// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "array_with_unit.h"
%}

%base_import(generic_object)
%import "unit.i"
%fp_shared_ptr(FullPhysics::ArrayWithUnit<double, 1>)
%fp_shared_ptr(FullPhysics::ArrayWithUnit<double, 2>)
%fp_shared_ptr(FullPhysics::ArrayWithUnit<double, 3>)

namespace FullPhysics {
template<class T, int D> class ArrayWithUnit : public GenericObject {
public:
  ArrayWithUnit();
  ArrayWithUnit(const ArrayWithUnit<T, D>& V);
  ArrayWithUnit(const blitz::Array<T, D>& FORCE_COPY, const Unit& Value_units);
  ArrayWithUnit(const blitz::Array<T, D>& FORCE_COPY, 
		const std::string& Value_units_name);
  %extend {
    // Need to return value through function to have 
    // conversion to numpy array happen
    blitz::Array<T, D> _value() const
    { return $self->value; }
    void _value(const blitz::Array<T, D>& Dat)
    { $self->value.reference(Dat.copy()); }
    Unit _units() const { return $self->units; }
    void _units(const Unit& U) { $self->units = U; }
  }

  ArrayWithUnit<T, D>& operator*=(const ArrayWithUnit<T, D>& V);
  ArrayWithUnit<T, D>& operator/=(const ArrayWithUnit<T, D>& V);
  ArrayWithUnit<T, D> convert(const Unit& R) const;
  ArrayWithUnit<T, D> convert(const std::string& R) const;
  ArrayWithUnit<T, D> convert_wave(const Unit& R) const;
  ArrayWithUnit<T, D> convert_wave(const std::string& R) const;
  %python_attribute(rows, int);
  %python_attribute(cols, int);
  %python_attribute(depth, int);
  std::string print_to_string() const;
  %pythoncode {
@property
def value(self):    
  return self._value()

@value.setter
def value(self,v):
  self._value(v)

@property
def units(self):    
  return self._units()

@units.setter
def units(self,v):
  self._units(v)
  }
};
}

%template(ArrayWithUnit_double_1) FullPhysics::ArrayWithUnit<double, 1>;
%template(ArrayWithUnit_double_2) FullPhysics::ArrayWithUnit<double, 2>;
%template(ArrayWithUnit_double_3) FullPhysics::ArrayWithUnit<double, 3>;
