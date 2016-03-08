#include "altitude.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Altitude)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<Altitude> >::*pbt1)(
        const std::vector<boost::shared_ptr<Altitude> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<Altitude> >, VectorAltitude)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<Altitude> >::push_back))
REGISTER_LUA_END()
#endif
