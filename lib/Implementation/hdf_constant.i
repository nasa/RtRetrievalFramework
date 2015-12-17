// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "hdf_constant.h"
%}

%base_import(constant)
%import "hdf_file.i"
%import "double_with_unit.i"
%fp_shared_ptr(FullPhysics::HdfConstant);

namespace FullPhysics {
class HdfConstant: public Constant {
public:
  HdfConstant(const boost::shared_ptr<HdfFile>& Hdf_file);
  %python_attribute_derived(rayleigh_depolarization_factor, double);
  %python_attribute_derived(rayleigh_a, DoubleWithUnit);
  %python_attribute_derived(rayleigh_b, DoubleWithUnit);
  %python_attribute_derived(molar_weight_dry_air, DoubleWithUnit);
  %python_attribute_derived(molar_weight_water, DoubleWithUnit);
  %python_attribute_derived(avogadro_constant, DoubleWithUnit);
};

}

