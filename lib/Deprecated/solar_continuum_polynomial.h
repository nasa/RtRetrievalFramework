#ifndef SOLAR_CONTINUUM_POLYNOMIAL_H
#define SOLAR_CONTINUUN_POLYNOMIAL_H
#include "solar_continuum_spectrum.h"
#include "array_with_unit.h"
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This class calculates the solar continuum spectrum.

  This particular implementation uses a polynomial parametrization
  to calculate the Solar Planck Function.
*******************************************************************/

class SolarContinuumPolynomial : public SolarContinuumSpectrum {
public:
  SolarContinuumPolynomial(const ArrayWithUnit<double, 1>& Param,
			   bool Convert_from_photon = true);
  virtual ~SolarContinuumPolynomial() {}
  virtual void print(std::ostream& Os) const;
  virtual Spectrum solar_continuum_spectrum(
     const SpectralDomain& spec_domain) const;
private:
  ArrayWithUnit<double, 1> param;
  bool convert_from_photon;
};
}
#endif
