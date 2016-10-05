#ifndef RADIANCE_SCALING_LINEAR_FIT_H
#define RADIANCE_SCALING_LINEAR_FIT_H
#include "radiance_scaling.h"
#include "double_with_unit.h"
#include "spectral_range.h"
#include "instrument_correction.h"

namespace FullPhysics {

/****************************************************************//**
  Implements a fitted radiance scaling correction where the correction
  is determined by a linear fit of the measured spectra to the radiance
  input instead of using the state vector for fitting.
*******************************************************************/

class RadianceScalingLinearFit : public RadianceScaling {
public:
//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Measured_radiance Measured radiance for the band 
/// \param Band_ref The wavenumber/wavelength that the scaling
///    coefficients are relative to.
/// \param Band_name Name of band
/// \param Do_offset If true then do the offset
//-----------------------------------------------------------------------
  
  RadianceScalingLinearFit(const SpectralRange& Measured_radiance,
                           const DoubleWithUnit& Band_ref,
                           const std::string& Band_name,
                           const bool Do_offset = true)
    : RadianceScaling(Band_ref, Band_name), measured_radiance(Measured_radiance),
      do_offset(Do_offset)
  { } 

  virtual ~RadianceScalingLinearFit() {}

  virtual boost::shared_ptr<InstrumentCorrection> clone() const;

  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const;

  virtual void print(std::ostream& Os) const;

  //-----------------------------------------------------------------------
  ///  
  //-----------------------------------------------------------------------

  blitz::Array<double, 1> radiance_scaling_coeff_uncertainty() const 
  { 
    blitz::Array<double, 1> res(scaling_coeff.rows());
    res = 0;
    return res;
  }
private:
  SpectralRange measured_radiance;
  bool do_offset;
};
}
#endif
