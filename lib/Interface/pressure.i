// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "pressure.h"
%}

%base_import(state_vector)
%import "array_ad_with_unit.i"
%import "auto_derivative_with_unit.i"

%fp_shared_ptr(FullPhysics::Pressure)
namespace FullPhysics {
  class Pressure;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Pressure>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Pressure>)

%template(ObservablePressure) FullPhysics::Observable<FullPhysics::Pressure>;
%template(ObserverPressure) FullPhysics::Observer<FullPhysics::Pressure>;

namespace FullPhysics {
// Allow this class to be derived from in Python.
%feature("director") Pressure;


class Pressure : virtual public StateVectorObserver, 
		 public Observable<Pressure> {
public:
  virtual ~Pressure();
  virtual void add_observer(Observer<Pressure>& Obs);
  virtual void remove_observer(Observer<Pressure>& Obs);
  %python_attribute(surface_pressure, AutoDerivativeWithUnit<double>)
  %python_attribute(surface_pressure_value, double)
  %python_attribute_abstract(pressure_grid, ArrayAdWithUnit<double, 1>)
  %python_attribute(number_layer, int)
  %python_attribute(number_level, int)
  %python_attribute(max_number_level, virtual int)
  virtual boost::shared_ptr<Pressure> clone() const = 0;
  std::string print_to_string() const;
};
}

