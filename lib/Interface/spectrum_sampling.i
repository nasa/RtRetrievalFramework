// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "spectrum_sampling.h"
%}
%base_import(generic_object)
%import "spectral_domain.i"
%import "spectrum.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::SpectrumSampling);

namespace FullPhysics {
class SpectrumSampling : public GenericObject {
public:
  virtual ~SpectrumSampling();
  std::string print_to_string() const;
  %python_attribute(number_spectrometer, int);
  virtual SpectralDomain spectral_domain(int spec_index,
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Ils_half_width) const = 0;
  virtual SpectralDomain spectral_domain_interpolated(int Spec_index, 
		 const SpectralDomain& Lowres_grid, 
	         const DoubleWithUnit& Ils_half_width) const;
  virtual bool need_interpolation(int Spec_index) const;
};
}
