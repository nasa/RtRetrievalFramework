#ifndef GROUND_LAMBERTIAN_OUTPUT_H
#define GROUND_LAMBERTIAN_OUTPUT_H
#include "register_output_base.h"
#include "ground_lambertian.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the GroundLambertian class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the GroundLambertian class.
*******************************************************************/
class GroundLambertianOutput : public RegisterOutputBase {
public:
  GroundLambertianOutput(const boost::shared_ptr<GroundLambertian>& Lambertian, const std::vector<std::string>& Hdf_band_names) 
    : lambertian(Lambertian), hdf_band_names(Hdf_band_names), surface_type("Lambertian") {}
  virtual ~GroundLambertianOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<GroundLambertian> lambertian;
  std::vector<std::string> hdf_band_names;
  std::string surface_type;
};

}
#endif
