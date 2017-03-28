#ifndef GROUND_BRDF_OUTPUT_H
#define GROUND_BRDF_OUTPUT_H
#include "register_output_base.h"
#include "ground_brdf.h"
#include "level_1b.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the GroundBrdf class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the GroundBrdf class.
*******************************************************************/
class GroundBrdfOutput : public RegisterOutputBase {
public:
  GroundBrdfOutput(const boost::shared_ptr<GroundBrdf>& Brdf, 
                   const boost::shared_ptr<Level1b>& L1b,
                   const std::vector<std::string>& Hdf_band_names) 
    : brdf(Brdf), l1b(L1b), hdf_band_names(Hdf_band_names) {}
  virtual ~GroundBrdfOutput() {}
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<GroundBrdf> brdf;
  boost::shared_ptr<Level1b> l1b;
  std::vector<std::string> hdf_band_names;
  mutable std::string surface_type;
};

}
#endif
