// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "constant.h"
%}
%base_import(generic_object)
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::Constant);

namespace FullPhysics {
class Constant : public GenericObject {
public:
  std::string print_to_string() const;
  %python_attribute_abstract(rayleigh_depolarization_factor, double);
  %python_attribute_abstract(rayleigh_a, DoubleWithUnit);
  %python_attribute_abstract(rayleigh_b, DoubleWithUnit);
  %python_attribute_abstract(molar_weight_dry_air, DoubleWithUnit);
  %python_attribute_abstract(molar_weight_water, DoubleWithUnit);
  %python_attribute_abstract(avogadro_constant, DoubleWithUnit);
};

}

