// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"
%{
#include "absco.h"
%}
%base_import(gas_absorption)
%import "double_with_unit.i"
%import "auto_derivative_with_unit.i"
%fp_shared_ptr(FullPhysics::Absco);
namespace FullPhysics {
class Absco : public GasAbsorption {
public:
  %python_attribute(number_broadener_vmr, int)
  %python_attribute(number_layer, int)
  %python_attribute(number_temperature, int)
  %python_attribute_abstract(broadener_vmr_grid,blitz::Array<double, 1>)
  %python_attribute_abstract(pressure_grid,blitz::Array<double, 1>)
  %python_attribute_abstract(temperature_grid,blitz::Array<double, 2>)
  %python_attribute(broadener_name, std::string)
  virtual double table_scale(double wn) const;
  virtual DoubleWithUnit absorption_cross_section(double Wn, 
     const DoubleWithUnit& Press, 
     const DoubleWithUnit& Temp,
     const DoubleWithUnit& Broadener_vmr) const;
  virtual AutoDerivativeWithUnit<double>
  absorption_cross_section(double wn, 
    const DoubleWithUnit& press, 
    const AutoDerivativeWithUnit<double>& temp,
    const AutoDerivativeWithUnit<double>& broadener_vmr) const;
  %extend {
    blitz::Array<double, 3> read_double(double wn) const 
    { return $self->read<double>(wn); }
    blitz::Array<float, 3> read_float(double wn) const 
    { return $self->read<float>(wn); }
  }
};
}

%template(vector_absco) std::vector<boost::shared_ptr<FullPhysics::Absco> >;
