// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "solar_absorption_gfit_file.h"
%}
%base_import(solar_absorption_spectrum)
%import "spectrum.i"
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::SolarAbsorptionGfitFile);

namespace FullPhysics {
class SolarAbsorptionGfitFile : public SolarAbsorptionSpectrum {
public:
  SolarAbsorptionGfitFile(const std::string& Line_list_file,
                          double Fraction_solar_diameter);
  virtual Spectrum solar_absorption_spectrum(
     const SpectralDomain& spec_domain) const;
  %python_attribute(fraction_solar_diameter, double)
};
}
