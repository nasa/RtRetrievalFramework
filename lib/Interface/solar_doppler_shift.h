#ifndef SOLAR_DOPPLER_SHIFT_H
#define SOLAR_DOPPLER_SHIFT_H
#include "printable.h"
#include "double_with_unit.h"
#include "spectral_domain.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This class handles the solar Doppler stretch to calculate the shift
  of the solar lines with respect to the telluric lines
*******************************************************************/
class SolarDopplerShift : public Printable<SolarDopplerShift> {
public:
  virtual ~SolarDopplerShift() {}

//-----------------------------------------------------------------------
/// Return Earth-Sun distance.
/// \return Earth-Sun distance.
//-----------------------------------------------------------------------
  
  virtual DoubleWithUnit solar_distance() const = 0;

//-----------------------------------------------------------------------
/// Shift wavenumbers to account for doppler stretch.
/// \param Spec_domain wavenumber/wavelength
/// \return Wavenumber/wavelength with Doppler stretch.
//-----------------------------------------------------------------------

  virtual SpectralDomain doppler_stretch(
     const SpectralDomain& Spec_domain) const = 0;
  virtual void print(std::ostream& Os) const {Os << "SolarDopplerShift";}
};
}
#endif
