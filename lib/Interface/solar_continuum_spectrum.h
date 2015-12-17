#ifndef SOLAR_CONTINUUM_SPECTRUM_H
#define SOLAR_CONTINUUM_SPECTRUM_H
#include "printable.h"
#include "spectrum.h"

namespace FullPhysics {
/****************************************************************//**
  This class calculates the solar continuum spectrum
*******************************************************************/
class SolarContinuumSpectrum : public Printable <SolarContinuumSpectrum> {
public:
  virtual ~SolarContinuumSpectrum() {}

//-----------------------------------------------------------------------
/// This calculate the solar continuum spectrum. This calculates this
/// at an Earth-Sun distance of 1 AU, this needs to be scaled by the
/// square of the actual distance.
///
/// \param Spec_domain Wavenumber/Wavelength to return solar continuum spectrum
///                    for. 
/// \return The solar continuum spectrum at 1 AU.
///
/// Note that the spectral domain here is in the solar rest frame,
/// *not* the earth rest frame used in most other places. The class
/// SolarAbsorptionAndContinuum handles this conversion internally.
//-----------------------------------------------------------------------

  virtual Spectrum solar_continuum_spectrum(
     const SpectralDomain& Spec_domain) const = 0;
  virtual void print(std::ostream& Os) const {Os << "SolarContinuumSpectrum";}
};
}
#endif
