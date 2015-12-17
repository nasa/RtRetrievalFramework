#include "forward_model_cost_function.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ForwardModelCostFunction, CostFunction)
.def(luabind::constructor<const boost::shared_ptr<ForwardModel>&>())
REGISTER_LUA_END()
#endif
