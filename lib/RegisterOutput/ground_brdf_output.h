#ifndef GROUND_BRDF_OUTPUT_H
#define GROUND_BRDF_OUTPUT_H
#include "register_output_base.h"
#include "ground_brdf.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the GroundBrdf class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the GroundBrdf class.
*******************************************************************/
class GroundBrdfOutput : public RegisterOutputBase {
public:
  GroundBrdfOutput(const boost::shared_ptr<GroundBrdf>& Brdf, const blitz::Array<double, 1>& Sza, const blitz::Array<double, 1>& Eff_alb_table_cos_sza, const blitz::Array<double, 1>& Eff_alb_table_intensity_scaling, const std::vector<std::string>& Hdf_band_names) 
    : brdf(Brdf), sza(Sza), eff_alb_table_cos_sza(Eff_alb_table_cos_sza), eff_alb_table_intensity_scaling(Eff_alb_table_intensity_scaling), hdf_band_names(Hdf_band_names) {}
  virtual ~GroundBrdfOutput() {}
  double interpolated_intensity_scaling(int spec_idx) const;
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<GroundBrdf> brdf;
  const blitz::Array<double, 1>& sza;
  const blitz::Array<double, 1>& eff_alb_table_cos_sza;
  const blitz::Array<double, 1>& eff_alb_table_intensity_scaling;
  std::vector<std::string> hdf_band_names;
  mutable std::string surface_type;
};

}
#endif
