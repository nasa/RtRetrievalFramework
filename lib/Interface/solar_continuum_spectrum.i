// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "solar_continuum_spectrum.h"
%}

%base_import(generic_object)
%import "spectral_domain.i"
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::SolarContinuumSpectrum);

namespace FullPhysics {
class SolarContinuumSpectrum : public GenericObject {
public:
  virtual ~SolarContinuumSpectrum();
  std::string print_to_string() const;
  virtual Spectrum solar_continuum_spectrum(
     const SpectralDomain& spec_domain) const = 0;
};
}
