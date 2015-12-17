// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "solar_absorption_table.h"
%}
%base_import(solar_absorption_spectrum)
%import "hdf_file.i"
%import "spectrum.i"
%import "spectral_domain.i"

%fp_shared_ptr(FullPhysics::SolarAbsorptionTable);

namespace FullPhysics {
class SolarAbsorptionTable : public SolarAbsorptionSpectrum {
public:
  SolarAbsorptionTable(const HdfFile& Hdf_static_input,
			 const std::string& Hdf_group);
  virtual Spectrum solar_absorption_spectrum(
     const SpectralDomain& spec_domain) const;
};
}
