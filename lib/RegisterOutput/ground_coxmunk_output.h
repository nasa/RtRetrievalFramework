#ifndef GROUND_COXMUNK_OUTPUT_H
#define GROUND_COXMUNK_OUTPUT_H
#include "register_output_base.h"
#include "ground_coxmunk.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the GroundCoxmunk class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the GroundCoxmunk class.
*******************************************************************/
class GroundCoxmunkOutput : public RegisterOutputBase {
public:
  GroundCoxmunkOutput(const boost::shared_ptr<GroundCoxmunk>& Coxmunk) 
    : coxmunk(Coxmunk), surface_type("Coxmunk") {}
  virtual ~GroundCoxmunkOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<GroundCoxmunk> coxmunk;
  std::string surface_type;
};
}
#endif
