#ifndef SOLAR_ABSORPTION_AND_CONTINUUM_H
#define SOLAR_ABSORPTION_AND_CONTINUUM_H
#include "solar_model.h"
#include "solar_doppler_shift.h"
#include "solar_absorption_spectrum.h"
#include "solar_continuum_spectrum.h"
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/****************************************************************//**
  This applies a solar model to radiances to model the incoming solar
  irradiance.

  This implementation is a common division of the solar model into

  -# A Doppler correction
  -# A solar absorption spectrum
  -# A solar continuum spectrum

  This uses 3 objects to do the work, a SolarDopplerShift, a
  SolarAbsorptionSpectrum, and a SolarContinuumSpectrum object. This
  class stitches these objects together to create the full spectrum.
*******************************************************************/

class SolarAbsorptionAndContinuum : public SolarModel {
public:
//-----------------------------------------------------------------------
/// Create a SolarModel that uses the given doppler shift, absorption
/// spectrum, and continuum spectrum.
//-----------------------------------------------------------------------

  SolarAbsorptionAndContinuum(
     const boost::shared_ptr<SolarDopplerShift>& doppler_shiftv, 
     const boost::shared_ptr<SolarAbsorptionSpectrum>& absorption_spectrumv,
     const boost::shared_ptr<SolarContinuumSpectrum>& continuum_spectrumv)
    : doppler_shift_(doppler_shiftv),
      absorption_spectrum_(absorption_spectrumv),
      continuum_spectrum_(continuum_spectrumv)
  {}

  virtual ~SolarAbsorptionAndContinuum() {}


  //-----------------------------------------------------------------------
  /// Clone a SolarAbsorptionAndContinuum object
  //----------------------------------------------------------------------- 

  virtual boost::shared_ptr<SpectrumEffect> clone() const {
    boost::shared_ptr<SpectrumEffect> res
      (new SolarAbsorptionAndContinuum(doppler_shift_,
				       absorption_spectrum_,
				       continuum_spectrum_));
    return res;
  }

//-----------------------------------------------------------------------
/// The SolarDopplerShift object used by this class.
//-----------------------------------------------------------------------

  const SolarDopplerShift& doppler_shift() const {return *doppler_shift_;}

//-----------------------------------------------------------------------
/// The SolarDopplerShift object used by this class, as a ptr.
//-----------------------------------------------------------------------

  const boost::shared_ptr<SolarDopplerShift>& doppler_shift_ptr() const 
  {return doppler_shift_;}

//-----------------------------------------------------------------------
/// The SolarAbsorptionSpectrum object used by this class.
//-----------------------------------------------------------------------

  const SolarAbsorptionSpectrum& absorption_spectrum() const 
  {return *absorption_spectrum_;}

//-----------------------------------------------------------------------
/// The SolarAbsorptionSpectrum object used by this class, as a ptr.
//-----------------------------------------------------------------------

  const boost::shared_ptr<SolarAbsorptionSpectrum>& absorption_spectrum_ptr() 
    const {return absorption_spectrum_;}

//-----------------------------------------------------------------------
/// The SolarContinuumSpectrum object used by this class.
//-----------------------------------------------------------------------

  const SolarContinuumSpectrum& continuum_spectrum() const 
  {return *continuum_spectrum_;}

//-----------------------------------------------------------------------
/// The SolarContinuumSpectrum object used by this class, as a ptr.
//-----------------------------------------------------------------------

  const boost::shared_ptr<SolarContinuumSpectrum>& continuum_spectrum_ptr() 
    const {return continuum_spectrum_;}
  virtual void print(std::ostream& Os) const;
  virtual Spectrum solar_spectrum(const SpectralDomain& Spec_domain) const;
private:
  boost::shared_ptr<SolarDopplerShift> doppler_shift_;
  boost::shared_ptr<SolarAbsorptionSpectrum> absorption_spectrum_;
  boost::shared_ptr<SolarContinuumSpectrum> continuum_spectrum_;
};
}
#endif
