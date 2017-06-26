#ifndef SOLAR_DOPPLER_SHIFT_POLYNOMIAL_H
#define SOLAR_DOPPLER_SHIFT_POLYNOMIAL_H
#include "solar_doppler_shift.h"
#include "fp_time.h"
#include "double_with_unit.h"
#include <boost/array.hpp>
#include "default_constant.h"

namespace FullPhysics {
/****************************************************************//**
  This class handles the solar Doppler stretch to calculate the shift
  of the solar lines with respect to the telluric lines.

  This implementation uses a fit based on a 6th order polynomial fit
  to data from http://eclipse.gsfc.nasa.gov/TYPE/TYPE.html
*******************************************************************/
class SolarDopplerShiftPolynomial : public SolarDopplerShift {
public:
  SolarDopplerShiftPolynomial(const Time& t, 
                              const DoubleWithUnit& lat, 
                              const DoubleWithUnit& sol_zen,
                              const DoubleWithUnit& sol_az, 
                              const DoubleWithUnit& elevation,
                              const Constant& constant = DefaultConstant(),
                              bool apply_doppler_shift = true);
  SolarDopplerShiftPolynomial(double Doppler_shift,
                              const Time& t, 
                              const Constant& constant = DefaultConstant(),
                              bool apply_doppler_shift = true);

  virtual ~SolarDopplerShiftPolynomial() {}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const;
  virtual DoubleWithUnit solar_distance() const {return solar_distance_;}
  virtual SpectralDomain doppler_stretch(
     const SpectralDomain& Spec_domain) const;

//-----------------------------------------------------------------------
/// Velocity of the center of the earth to the center of the sun. 
/// Positive means they are getting farther apart. This does *not*
/// include the rotation of the earth.
//-----------------------------------------------------------------------

  DoubleWithUnit solar_velocity() const {return solar_velocity_;}

//-----------------------------------------------------------------------
/// Total velocity, including rotation of earth.
//-----------------------------------------------------------------------

  DoubleWithUnit total_velocity() const
  { return solar_velocity_ + doppler_rot_earth_sun_; }

  double doppler_shift() const {return doppler_shift_;}
  
private:
  DoubleWithUnit solar_distance_; ///< Earth-Sun distance
  DoubleWithUnit solar_velocity_;
  DoubleWithUnit doppler_rot_earth_sun_; 
  double doppler_shift_;         ///< Doppler shift factor.
  bool apply_doppler_shift_;         ///< If true, apply doppler shift.
  bool calculated_doppler_shift;
  void calc_solar_distance(const Constant& constant, const Time& t);        
};

}
#endif
