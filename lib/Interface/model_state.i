// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "model_state.h"
%}
%base_import(problem_state)
%fp_shared_ptr(FullPhysics::ModelState);

namespace FullPhysics {
class ModelState : virtual public ProblemState {
public:
  ModelState();
  ModelState(const ModelState& s);
  virtual ~ModelState();
  virtual void set(const ModelState& s);
  virtual void clear();
};
}
