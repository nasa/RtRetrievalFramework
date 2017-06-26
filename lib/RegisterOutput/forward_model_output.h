#ifndef FORWARD_MODEL_OUTPUT_H
#define FORWARD_MODEL_OUTPUT_H
#include "register_output_base.h"
#include "forward_model.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the ForwardModel class that should be
  written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the ForwardModel class.
*******************************************************************/
class ForwardModelOutput : public RegisterOutputBase {
public:
    ForwardModelOutput(const boost::shared_ptr<ForwardModel>&
                       Fm)
        : fm(Fm) {}
    virtual ~ForwardModelOutput() {}
    virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
    boost::shared_ptr<ForwardModel> fm;
};
}
#endif
