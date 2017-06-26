// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "spectrum_sampling_fixed_spacing.h"
%}
%base_import(spectrum_sampling)
%import "array_with_unit.i"
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::SpectrumSamplingFixedSpacing);
namespace FullPhysics {
class SpectrumSamplingFixedSpacing : public SpectrumSampling {
public:
  SpectrumSamplingFixedSpacing(const ArrayWithUnit<double, 1>& Spec_spacing);
  virtual SpectralDomain spectral_domain(int spec_index,
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Ils_half_width) const;
};
}
