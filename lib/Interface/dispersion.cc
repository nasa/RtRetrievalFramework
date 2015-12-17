#include "dispersion.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_CLASS(Dispersion)
.def("pixel_grid", &Dispersion::pixel_grid)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<Dispersion> >,
			VectorDispersion)
.def(luabind::constructor<>())
.def("push_back", &std::vector<boost::shared_ptr<Dispersion> >::push_back)
REGISTER_LUA_END()

#endif
