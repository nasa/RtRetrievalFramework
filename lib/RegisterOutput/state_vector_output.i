// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "state_vector_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "state_vector.i"

%fp_shared_ptr(FullPhysics::StateVectorOutput);

namespace FullPhysics {
class StateVectorOutput : public RegisterOutputBase {
public:
  StateVectorOutput(const boost::shared_ptr<StateVector>& Sv);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
};
}


