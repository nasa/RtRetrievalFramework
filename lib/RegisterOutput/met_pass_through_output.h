#ifndef MET_PASS_THROUGH_H
#define MET_PASS_THROUGH_H

#include "register_output_base.h"
#include "meteorology.h"

namespace FullPhysics {

/****************************************************************//**
  Writes source filenames into the output file
*******************************************************************/
class MetPassThroughOutput : public RegisterOutputBase {
public:
    MetPassThroughOutput(const boost::shared_ptr<Meteorology>& met) : met_(met) {}
    virtual ~MetPassThroughOutput() {}
    virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
    const boost::shared_ptr<Meteorology> met_;
};
}
#endif
