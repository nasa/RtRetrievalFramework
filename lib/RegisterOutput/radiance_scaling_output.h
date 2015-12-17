#ifndef RADIANCE_SCALING_OUTPUT_H
#define RADIANCE_SCALING_OUTPUT_H
#include "register_output_base.h"
#include "radiance_scaling.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the RadianceScaling class that
  should be written as output.
*******************************************************************/
class RadianceScalingOutput : public RegisterOutputBase {
public:
  RadianceScalingOutput(const boost::shared_ptr<RadianceScaling>& Rad_scaling,
                        const std::string& Hdf_band_name) 
    : rad_scaling(Rad_scaling),
      hdf_band_name(Hdf_band_name) {}
  virtual ~RadianceScalingOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<RadianceScaling> rad_scaling;
  std::string hdf_band_name;
};
}
#endif
