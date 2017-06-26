// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "radiance_scaling_output.h"
#include "sub_state_vector_array.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "radiance_scaling.i"

%fp_shared_ptr(FullPhysics::RadianceScalingOutput);

namespace FullPhysics {
class RadianceScalingOutput : public RegisterOutputBase {
public:
  RadianceScalingOutput(const boost::shared_ptr<RadianceScaling>& Rad_scaling,
                        const std::string& Hdf_band_name);
  virtual ~RadianceScalingOutput();
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}
