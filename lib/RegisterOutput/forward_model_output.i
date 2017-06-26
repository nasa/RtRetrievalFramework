// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "forward_model_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "forward_model.i"

%fp_shared_ptr(FullPhysics::ForwardModelOutput);

namespace FullPhysics {
class ForwardModelOutput : public RegisterOutputBase {
public:
  ForwardModelOutput(const boost::shared_ptr<ForwardModel>& Fm);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
};
}


