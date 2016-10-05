// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "default_constant.h"
%}

%base_import(constant);

%fp_shared_ptr(FullPhysics::DefaultConstant);

namespace FullPhysics {
class DefaultConstant: public Constant {
public:
  DefaultConstant();
  %python_attribute(rayleigh_depolarization_factor, virtual double);
  %python_attribute(rayleigh_a, virtual DoubleWithUnit);
  %python_attribute(rayleigh_b, virtual DoubleWithUnit);
  %python_attribute(molar_weight_dry_air, virtual DoubleWithUnit);
  %python_attribute(molar_weight_water, virtual DoubleWithUnit);
  %python_attribute(avogadro_constant, virtual DoubleWithUnit);
};

}

