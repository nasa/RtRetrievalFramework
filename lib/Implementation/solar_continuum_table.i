// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "solar_continuum_table.h"
%}
%base_import(solar_continuum_spectrum)
%import "hdf_file.i"
%import "spectrum.i"
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::SolarContinuumTable);

namespace FullPhysics {
class SolarContinuumTable : public SolarContinuumSpectrum {
public:
  SolarContinuumTable(const HdfFile& F, const std::string& Hdf_group,
		      bool Convert_from_photon = true);
  virtual Spectrum solar_continuum_spectrum(
     const SpectralDomain& spec_domain) const;
};
}
