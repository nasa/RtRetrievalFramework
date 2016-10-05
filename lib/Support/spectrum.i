// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"

%{
#include "spectrum.h"
%}
%base_import(generic_object)
%import "spectral_domain.i"
%import "spectral_range.i"

%fp_shared_ptr(FullPhysics::Spectrum)
namespace FullPhysics {
class Spectrum : public GenericObject {
public:
  std::string print_to_string() const;
  Spectrum(const SpectralDomain& Spec_domain, 
	   const SpectralRange& Spec_range);
  %python_attribute(spectral_domain, SpectralDomain)
  %python_attribute(spectral_range, SpectralRange)
  %pythoncode {
@property
def wavenumber(self):
    return self.spectral_domain.wavenumber()

@property
def wavelength(self):
    return self.spectral_domain.wavelength()

@property
def value(self):
    return self.spectral_range.data

@property
def units(self):
    return self.spectral_range.units

}
};
}

%template(vector_spectrum) std::vector<FullPhysics::Spectrum>;
