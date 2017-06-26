// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "solar_doppler_shift_l1b.h"
%}
%base_import(solar_doppler_shift)
%import "double_with_unit.i"
%import "spectral_domain.i"

%fp_shared_ptr(FullPhysics::SolarDopplerShiftL1b);

namespace FullPhysics {
class SolarDopplerShiftL1b : public SolarDopplerShift {
public:
  SolarDopplerShiftL1b(const DoubleWithUnit& Solar_distance, 
		       const DoubleWithUnit& Solar_relative_velocity,
		       bool Apply_doppler_shift = true);
  %python_attribute(doppler_shift, double)
  %python_attribute(solar_relative_velocity, DoubleWithUnit)
  %python_attribute(solar_distance, DoubleWithUnit)
  virtual SpectralDomain doppler_stretch(
     const SpectralDomain& Spec_domain) const;
};
}
