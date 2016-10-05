#include <model_state.h>


using namespace FullPhysics;


void ModelState::set(const ModelState& s)
{
  ProblemState::set(s);
  M.reference(s.M.copy());
  K.reference(s.K.copy());
}
