#ifndef ECMWF_PASS_THROUGH_H
#define ECMWF_PASS_THROUGH_H

#include "register_output_base.h"
#include "ecmwf.h"

namespace FullPhysics {

/****************************************************************//**
  Writes source filenames into the output file
*******************************************************************/
class EcmwfPassThroughOutput : public RegisterOutputBase {
public:
    EcmwfPassThroughOutput(const boost::shared_ptr<Ecmwf>& ecmwf) : ecmwf_(ecmwf) {}
    virtual ~EcmwfPassThroughOutput() {}
    virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
    const boost::shared_ptr<Ecmwf> ecmwf_;
};
}
#endif
