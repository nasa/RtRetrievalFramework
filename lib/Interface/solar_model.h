#ifndef SOLAR_MODEL_H
#define SOLAR_MODEL_H
#include "printable.h"
#include "array_ad.h"
#include "spectrum.h"
#include "spectrum_effect.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This applies a solar model to reflectance to model the incoming solar
  irradiance.
*******************************************************************/

class SolarModel : public SpectrumEffect {
public:
  virtual ~SolarModel() {}
  virtual Spectrum apply_solar_model(const Spectrum& Spec) const;

//-----------------------------------------------------------------------
/// Calculate solar spectrum.
///
/// \param Spec_domain Wavenumber/Wavelength reflectance is given
/// \return Solar spectrum. This should have units commensurate with
/// something like W / m^2 / cm^-1.
///
/// Note that the wavenumber/frequency are in the earth rest
/// frame. The solar model may need to work in the solar rest frame,
/// bu the conversion to this is internal. The input and output from
/// this function should be in the earth rest frame.
//-----------------------------------------------------------------------

  virtual Spectrum solar_spectrum(const SpectralDomain& Spec_domain) const = 0;
  virtual void print(std::ostream& Os) const {Os << "SolarModel";}

  virtual void apply_effect(Spectrum& Spec, 
		    const ForwardModelSpectralGrid& Forward_model_grid) const {
    Spec = apply_solar_model(Spec);
  }

  virtual std::string name() const { return "solar_model"; }
};
}
#endif
