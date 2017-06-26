// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "oco_forward_model_output.h"
#include "ils_instrument.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "oco_forward_model.i"

%fp_shared_ptr(FullPhysics::OcoForwardModelOutput);

namespace FullPhysics {
class OcoForwardModelOutput : public RegisterOutputBase {
public:
  OcoForwardModelOutput(const boost::shared_ptr<OcoForwardModel>& Fm);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
};
}
