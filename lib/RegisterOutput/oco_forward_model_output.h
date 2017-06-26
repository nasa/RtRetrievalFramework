#ifndef OCO_FORWARD_MODEL_OUTPUT_H
#define OCO_FORWARD_MODEL_OUTPUT_H
#include "register_output_base.h"
#include "oco_forward_model.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the OcoForwardModel class that 
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just 
  part of the OcoForwardModel class.
*******************************************************************/
class OcoForwardModelOutput : public RegisterOutputBase {
public:
    OcoForwardModelOutput(const boost::shared_ptr<OcoForwardModel>&
                       Fm)
        : fm(Fm) {}
    virtual ~OcoForwardModelOutput() {}
    virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
    boost::shared_ptr<OcoForwardModel> fm;
};
}
#endif
