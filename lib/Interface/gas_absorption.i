// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include <std_vector.i>
%include "common.i"

%{
#include "gas_absorption.h"
%}
%base_import(generic_object)

%import "double_with_unit.i"
%import "auto_derivative_with_unit.i"

%fp_shared_ptr(FullPhysics::GasAbsorption);
namespace FullPhysics {
class GasAbsorption : public GenericObject {
public:
  virtual ~GasAbsorption();
  std::string print_to_string() const;
  virtual bool have_data(double wn) const = 0;
  %python_attribute(broadener_name, virtual std::string);
  virtual DoubleWithUnit absorption_cross_section(double Wn, 
     const DoubleWithUnit& Press, 
     const DoubleWithUnit& Temp,
     const DoubleWithUnit& Broadener_vmr) const = 0;
  virtual AutoDerivativeWithUnit<double>
  absorption_cross_section(double wn, 
    const DoubleWithUnit& press, 
    const AutoDerivativeWithUnit<double>& temp,
    const AutoDerivativeWithUnit<double>& broadener_vmr) const = 0;
};
}

%template(vector_gas_absorption) std::vector<boost::shared_ptr<FullPhysics::GasAbsorption> >;
