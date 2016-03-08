#include "register_output_base.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(RegisterOutputBase)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<RegisterOutputBase> >::*pbt1)(
        const std::vector<boost::shared_ptr<RegisterOutputBase> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<RegisterOutputBase> >, VectorRegisterOutput)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<RegisterOutputBase> >::push_back))
REGISTER_LUA_END()
#endif
