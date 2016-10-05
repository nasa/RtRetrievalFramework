// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "spectral_window_range.h"
%}
%base_import(spectral_window)
%import "array_with_unit.i"
%import "heritage_file.i"
%import "double_with_unit.i"
%fp_shared_ptr(FullPhysics::SpectralWindowRange);

namespace FullPhysics {
class SpectralWindowRange : public SpectralWindow {
public:
  SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges);
  virtual std::vector<int> grid_indexes(const SpectralDomain& Grid, int Spec_index) const;
  %python_attribute(number_spectrometer, int)
  %python_attribute_with_set(range_array, ArrayWithUnit<double, 3>)
  %python_attribute_with_set(bad_sample_mask, blitz::Array<bool, 2>)
};
}
