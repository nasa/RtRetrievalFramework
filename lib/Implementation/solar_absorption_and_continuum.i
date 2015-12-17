// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "solar_absorption_and_continuum.h"
#include "sub_state_vector_array.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}
%base_import(solar_model)
%import "solar_doppler_shift.i"
%import "solar_absorption_spectrum.i"
%import "solar_continuum_spectrum.i"
%import "spectrum.i"
%import "spectral_domain.i"

%fp_shared_ptr(FullPhysics::SolarAbsorptionAndContinuum);

namespace FullPhysics {
class SolarAbsorptionAndContinuum : public SolarModel {
public:
  virtual ~SolarAbsorptionAndContinuum();
  SolarAbsorptionAndContinuum(
     const boost::shared_ptr<SolarDopplerShift>& doppler_shiftv, 
     const boost::shared_ptr<SolarAbsorptionSpectrum>& absorption_spectrumv,
     const boost::shared_ptr<SolarContinuumSpectrum>& continuum_spectrumv);
  virtual boost::shared_ptr<SpectrumEffect> clone() const;
  %python_attribute2(doppler_shift, doppler_shift_ptr, 
		     boost::shared_ptr<SolarDopplerShift>)
  %python_attribute2(absorption_spectrum, absorption_spectrum_ptr, 
		     boost::shared_ptr<SolarAbsorptionSpectrum>)
  %python_attribute2(continuum_spectrum, continuum_spectrum_ptr, 
		     boost::shared_ptr<SolarContinuumSpectrum>)
  virtual Spectrum solar_spectrum(const SpectralDomain& Spec_domain) const;
};
}
