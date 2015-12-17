// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "spectral_bound.h"
%}

%base_import(generic_object)
%import "array_with_unit.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::SpectralBound);

namespace FullPhysics {
class SpectralBound  : public GenericObject {
public:
  SpectralBound(const std::vector<DoubleWithUnit>& Lower_bound,
		const std::vector<DoubleWithUnit>& Upper_bound);
  SpectralBound(const ArrayWithUnit<double, 2>& Bound);
  std::string print_to_string() const;
  %python_attribute(number_spectrometer, int);
  DoubleWithUnit center(int Spec_index, const Unit& U) const;
  DoubleWithUnit lower_bound(int Spec_index) const;
  DoubleWithUnit upper_bound(int Spec_index) const;
  DoubleWithUnit lower_bound(int Spec_index, const Unit& U) const;
  DoubleWithUnit upper_bound(int Spec_index, const Unit& U) const;
  int spectral_index(const DoubleWithUnit& W) const;
};
}
