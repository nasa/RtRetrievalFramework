// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "solar_continuum_polynomial.h"
%}

%base_import(solar_continuum_spectrum)

%fp_shared_ptr(FullPhysics::SolarContinuumPolynomial);

namespace FullPhysics {
class SolarContinuumPolynomial : public SolarContinuumSpectrum {
public:
  SolarContinuumPolynomial(const ArrayWithUnit<double, 1>& Param,
			   bool Convert_from_photon = true);
  virtual Spectrum solar_continuum_spectrum(
     const SpectralDomain& spec_domain) const;
};
}
