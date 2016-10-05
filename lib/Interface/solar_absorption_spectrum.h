#ifndef SOLAR_ABSORPTION_SPECTRUM_H
#define SOLAR_ABSORPTION_SPECTRUM_H
#include "printable.h"
#include "spectrum.h"

namespace FullPhysics {
/****************************************************************//**
  This class calculates the solar absorption spectrum.
*******************************************************************/

class SolarAbsorptionSpectrum : public Printable<SolarAbsorptionSpectrum> {
public:
  virtual ~SolarAbsorptionSpectrum() {}

//-----------------------------------------------------------------------
/// This calculates the solar absorption spectrum.
///
/// \param Spec_domain Wavenumber/Wavelength to return solar
///        absorption spectrum for. 
/// \return The solar absorption spectrum. This is unit less, it is a
///         scaled value with 1.0 being no absorption.
///
/// Note that the spectral domain here is in the solar rest frame,
/// *not* the earth rest frame used in most other places. The class
/// SolarAbsorptionAndContinuum handles this conversion internally.
//-----------------------------------------------------------------------
  virtual Spectrum solar_absorption_spectrum(
     const SpectralDomain& Spec_domain) const = 0;
  virtual void print(std::ostream& Os) const {Os << "SolarAbsorptionSpectrum";}
};
}

#endif
