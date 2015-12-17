// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "solar_absorption_spectrum.h"
%}

%base_import(generic_object)
%import "spectral_domain.i"
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::SolarAbsorptionSpectrum);

namespace FullPhysics {
class SolarAbsorptionSpectrum : public GenericObject {
public:
  virtual ~SolarAbsorptionSpectrum();
  std::string print_to_string() const;
  virtual Spectrum solar_absorption_spectrum(
     const SpectralDomain& spec_domain) const = 0;
};
}
