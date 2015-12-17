#ifndef SOLAR_DOPPLER_SHIFT_L1B_H
#define SOLAR_DOPPLER_SHIFT_L1B_H
#include "solar_doppler_shift.h"
#include "fp_time.h"
#include "double_with_unit.h"
#include <boost/array.hpp>
#include "default_constant.h"

namespace FullPhysics {
/****************************************************************//**
  This class handles the solar Doppler stretch to calculate the shift
  of the solar lines with respect to the telluric lines.

  This implementation gets the Level 1b solar velocity and solar
  distance passed to it.
*******************************************************************/
class SolarDopplerShiftL1b : public SolarDopplerShift {
public:
  SolarDopplerShiftL1b(const DoubleWithUnit& Solar_distance, 
		       const DoubleWithUnit& Solar_relative_velocity,
		       bool Apply_doppler_shift = true);

  virtual ~SolarDopplerShiftL1b() {}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const;
  virtual DoubleWithUnit solar_distance() const {return solar_distance_;}
  virtual SpectralDomain doppler_stretch(
     const SpectralDomain& Spec_domain) const;

//-----------------------------------------------------------------------
/// Velocity of the sounding to the center of the sun.  Positive means
/// they are getting farther apart. This does include the rotation of 
/// the earth (cf solar_velocity in SolarDopplerShiftPolynomial which
/// does *not* include rotation of earth)
//-----------------------------------------------------------------------

  DoubleWithUnit solar_relative_velocity() const 
  {return solar_relative_velocity_;}

  double doppler_shift() const {return doppler_shift_;}
  
private:
  DoubleWithUnit solar_distance_; ///< Earth-Sun distance
  DoubleWithUnit solar_relative_velocity_;
  double doppler_shift_;         ///< Doppler shift factor.
  bool apply_doppler_shift_;         ///< If true, apply doppler shift.
};

}
#endif
