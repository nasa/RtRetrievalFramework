// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "solar_absorption_oco_file.h"
%}
%base_import(solar_absorption_spectrum)
%import "hdf_file.i"
%fp_shared_ptr(FullPhysics::SolarAbsorptionOcoFile);

namespace FullPhysics {
class SolarAbsorptionOcoFile : public SolarAbsorptionSpectrum {
public:
  SolarAbsorptionOcoFile(const HdfFile& Hdf_static_input,
			 const std::string& Hdf_group,
			 double Fraction_solar_diameter);
  virtual Spectrum solar_absorption_spectrum(
     const SpectralDomain& spec_domain) const;
  %python_attribute(number_line, int)
  %python_attribute(fraction_solar_diameter, double)
};
}
