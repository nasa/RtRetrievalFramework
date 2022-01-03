#ifndef GROUND_COXMUNK_SCALED_OUTPUT_H
#define GROUND_COXMUNK_SCALED_OUTPUT_H
#include "register_output_base.h"
#include "level_1b.h"
#include "ground_coxmunk_output.h"
#include "ground_coxmunk_scaled.h"
#include "ground_brdf_weight_output.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the GroundCoxmunkScaled
  class that should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the GroundCoxmunkScaled class.
*******************************************************************/
class GroundCoxmunkScaledOutput : public RegisterOutputBase {
public:
  GroundCoxmunkScaledOutput(const boost::shared_ptr<Level1b>& L1b,
                            const boost::shared_ptr<GroundCoxmunk>& Coxmunk,
		            const boost::shared_ptr<GroundBrdfWeight>& Brdf_weight,
                            const std::vector<std::string>& Hdf_band_names);
  virtual ~GroundCoxmunkScaledOutput() {}

  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  std::vector<std::string> hdf_band_names;
  const std::string surface_type;
  boost::shared_ptr<Level1b> l1b;
  boost::shared_ptr<GroundCoxmunk> coxmunk;
  boost::shared_ptr<GroundBrdfWeight> brdf_weight;
  boost::shared_ptr<GroundCoxmunkOutput> coxmunk_output;
  boost::shared_ptr<GroundBrdfWeightOutput> brdf_weight_output;
};
}
#endif
