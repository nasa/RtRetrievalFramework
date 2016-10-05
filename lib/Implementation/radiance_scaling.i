// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "radiance_scaling.h"
#include "sub_state_vector_array.h"
%}
%base_import(instrument_correction)
%import "double_with_unit.i"
%import "auto_derivative.i"
%import "array_ad.i"
%import "spectral_domain.i"
%import "spectral_range.i"
%fp_shared_ptr(FullPhysics::RadianceScaling);

namespace FullPhysics {
class RadianceScaling : virtual public InstrumentCorrection {
public:
  RadianceScaling(const DoubleWithUnit& Band_ref,
                  const std::string& Band_name);
  RadianceScaling(const ArrayAd<double, 1>& Scaling_coeff,
                  const DoubleWithUnit& Band_ref,
                  const std::string& Band_name);
  RadianceScaling(const ArrayAd<double, 1>& Scaling_coeff,
                  const AutoDerivative<double> Offset,
                  const DoubleWithUnit& Band_ref,
                  const std::string& Band_name); 
  virtual ~RadianceScaling();
  virtual void print(std::ostream& Os) const;
  virtual void apply_scaling(const SpectralDomain& Grid, SpectralRange& Radiance) const;
  %python_attribute(radiance_scaling_coeff, blitz::Array<double, 1>)
  %python_attribute(radiance_scaling_coeff_uncertainty, blitz::Array<double, 1>)
  %python_attribute(radiance_offset, double)
};
}
