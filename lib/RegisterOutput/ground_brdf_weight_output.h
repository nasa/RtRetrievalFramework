#ifndef GROUND_BRDF_WEIGHT_OUTPUT_H
#define GROUND_BRDF_WEIGHT_OUTPUT_H
#include "register_output_base.h"
#include "ground_brdf_weight.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the GroundBrdfWeight class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the GroundBrdfWeight class.
*******************************************************************/
class GroundBrdfWeightOutput : public RegisterOutputBase {
public:
    GroundBrdfWeightOutput(const boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, const std::vector<std::string>& Hdf_band_names)
      : brdf_weight(Brdf_weight), hdf_band_names(Hdf_band_names), surface_type("BRDFWeight") {}
    virtual ~GroundBrdfWeightOutput() {}

    double weight_intercept(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx);
    double weight_slope(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx);
    double weight_coeff(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx, int i);
    double weight_intercept_uncert(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx);
    double weight_slope_uncert(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx);
    double weight_coeff_uncert(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx, int i);

    virtual void register_output(const boost::shared_ptr<Output>& out) const;
    virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
    boost::shared_ptr<GroundBrdfWeight> brdf_weight;
    std::vector<std::string> hdf_band_names;
    std::string surface_type;
};

}
#endif
