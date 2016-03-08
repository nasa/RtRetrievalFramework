#include "gas_absorption.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(GasAbsorption)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<GasAbsorption> >::*pbt1)(
        const std::vector<boost::shared_ptr<GasAbsorption> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<GasAbsorption> >, VectorGasAbsorption)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<GasAbsorption> >::push_back))
REGISTER_LUA_END()
#endif
