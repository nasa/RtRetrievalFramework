#ifndef GROUND_LAMBERTIAN_PLUS_LAMB_OUTPUT_H
#define GROUND_LAMBERTIAN_PLUS_LAMB_OUTPUT_H
#include "register_output_base.h"
#include "ground_coxmunk_output.h"
#include "ground_lambertian_output.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the GroundCoxmunkPlusLambertian
  class that should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the GroundCoxmunkPlusLambertian class.
*******************************************************************/
class GroundCoxmunkPlusLambertianOutput : public RegisterOutputBase {
public:
  GroundCoxmunkPlusLambertianOutput(const boost::shared_ptr<GroundCoxmunk>& Coxmunk,
		                    const boost::shared_ptr<GroundLambertian>& Lambertian,
                                    const std::vector<std::string>& Hdf_band_names);
  virtual ~GroundCoxmunkPlusLambertianOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  const std::string surface_type;
  boost::shared_ptr<GroundCoxmunkOutput> coxmunk_output;
  boost::shared_ptr<GroundLambertianOutput> lambertian_output;
};
}
#endif
