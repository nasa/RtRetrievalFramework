// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "radiance_scaling_linear_fit.h"
#include "sub_state_vector_array.h"
%}
%base_import(radiance_scaling)
%import "double_with_unit.i"
%import "spectral_domain.i"
%import "spectral_range.i"
%fp_shared_ptr(FullPhysics::RadianceScalingLinearFit);

namespace FullPhysics {
class RadianceScalingLinearFit : public RadianceScaling {
public:
  RadianceScalingLinearFit(const SpectralRange& Measured_radiance,
                           const DoubleWithUnit& Band_ref,
                           const std::string& Band_name,
                           const bool Do_offset = true);
  virtual ~RadianceScalingLinearFit();
  virtual boost::shared_ptr<InstrumentCorrection> clone() const;
  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const;
  virtual void print(std::ostream& Os) const;
  %python_attribute(radiance_scaling_coeff_uncertainty, blitz::Array<double, 1>)
};
}
