// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "solar_doppler_shift_polynomial.h"
%}
%base_import(solar_doppler_shift)
%import "fp_time.i"
%import "double_with_unit.i"
%import "spectral_domain.i"
%import "constant.i"

%fp_shared_ptr(FullPhysics::SolarDopplerShiftPolynomial);

namespace FullPhysics {
class SolarDopplerShiftPolynomial : public SolarDopplerShift {
public:
  SolarDopplerShiftPolynomial(const Time& t, 
                              const DoubleWithUnit& lat, 
                              const DoubleWithUnit& sol_zen,
                              const DoubleWithUnit& sol_az, 
                              const DoubleWithUnit& elevation);
  SolarDopplerShiftPolynomial(double Doppler_shift,
                              const Time& t);
  SolarDopplerShiftPolynomial(const Time& t, 
                              const DoubleWithUnit& lat, 
                              const DoubleWithUnit& sol_zen,
                              const DoubleWithUnit& sol_az, 
                              const DoubleWithUnit& elevation,
                              const Constant& constant,
                              bool apply_doppler_shift = true);
  SolarDopplerShiftPolynomial(double Doppler_shift,
                              const Time& t, 
                              const Constant& constant,
                              bool apply_doppler_shift = true);
  %python_attribute(doppler_shift, double)
  %python_attribute(solar_velocity, DoubleWithUnit)
  %python_attribute(total_velocity, DoubleWithUnit)
  %python_attribute(solar_distance, DoubleWithUnit)
  virtual SpectralDomain doppler_stretch(
     const SpectralDomain& Spec_domain) const;
};
}
