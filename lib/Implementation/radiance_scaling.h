#ifndef RADIANCE_SCALING_H
#define RADIANCE_SCALING_H

#include "instrument_correction.h"
#include "double_with_unit.h"
#include "spectral_domain.h"
#include "spectral_range.h"
#include "printable.h"

namespace FullPhysics {

/****************************************************************//**
 This abstract class provides the generic capabilities for 
 applying a radiance scaling to a Radiance.

 The radiance scaling slope is referenced to a reference
 wavelength/wavenumber

 This class can support both a scale factor polynomial and a single
 radiance wide offset.
*******************************************************************/

class RadianceScaling : virtual public InstrumentCorrection,
                        public Printable<RadianceScaling> {

public:

  RadianceScaling(const DoubleWithUnit& Band_ref,
                  const std::string& Band_name) 
    : offset(0), band_ref(Band_ref), band_name(Band_name)
  {
    // Put in a sane default value
    scaling_coeff.resize(2, 0);
    scaling_coeff(0) = 1.0;
    scaling_coeff(1) = 0.0;
  }

  // Create a RadianceScaling
  RadianceScaling(const ArrayAd<double, 1>& Scaling_coeff,
                  const DoubleWithUnit& Band_ref,
                  const std::string& Band_name) 
    : scaling_coeff(Scaling_coeff),
      offset(0), band_ref(Band_ref), band_name(Band_name)
  { 
  }

  RadianceScaling(const ArrayAd<double, 1>& Scaling_coeff,
                  const AutoDerivative<double> Offset,
                  const DoubleWithUnit& Band_ref,
                  const std::string& Band_name) 
    : scaling_coeff(Scaling_coeff),
      offset(Offset), band_ref(Band_ref), band_name(Band_name)
  { 
  }

  virtual ~RadianceScaling() {};

  virtual void print(std::ostream& Os) const;

  /// Apply scaling and offset coefficients to Radiance 
  virtual void apply_scaling(const SpectralDomain& Grid, SpectralRange& Radiance) const;

  //-----------------------------------------------------------------------
  /// Return radiance scaling coefficients for the output file
  //-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> radiance_scaling_coeff() const { return scaling_coeff.value(); }

  //-----------------------------------------------------------------------
  /// Return radiance scaling coefficients uncertainty for the output file
  //-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> radiance_scaling_coeff_uncertainty() const = 0;

  //-----------------------------------------------------------------------
  /// Return radiance scaling offset for the output file
  //-----------------------------------------------------------------------

  virtual double radiance_offset() const { return offset.value(); }

protected:
  mutable ArrayAd<double, 1> scaling_coeff;
  mutable AutoDerivative<double> offset;
  DoubleWithUnit band_ref;
  std::string band_name;
};
}
#endif
