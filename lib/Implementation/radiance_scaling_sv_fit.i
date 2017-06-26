// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "radiance_scaling_sv_fit.h"
#include "sub_state_vector_array.h"
%}
%base_import(radiance_scaling)
%base_import(instrument_correction)
%base_import(sub_state_vector_array)
%import "double_with_unit.i"
%import "auto_derivative.i"
%import "array_ad.i"
%import "spectral_domain.i"
%import "spectral_range.i"

%fp_shared_ptr(FullPhysics::RadianceScalingSvFit);

namespace FullPhysics {
class RadianceScalingSvFit : public RadianceScaling, 
                              public SubStateVectorArray<InstrumentCorrection> {
public:
  RadianceScalingSvFit(const blitz::Array<double, 1>& Coeff, 
                        const blitz::Array<bool, 1>& Used_flag,
                        const DoubleWithUnit& Band_ref,
                        const std::string& Band_name);
  virtual ~RadianceScalingSvFit();
  virtual std::string state_vector_name_i(int i) const;
  virtual boost::shared_ptr<InstrumentCorrection> clone() const;
  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const;
  virtual void print(std::ostream& Os) const;
  virtual void notify_update(const StateVector& Sv);
  %python_attribute(radiance_scaling_coeff_uncertainty, blitz::Array<double, 1>)
};
}

