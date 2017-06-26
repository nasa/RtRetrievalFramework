#include "dispersion.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_CLASS(Dispersion)
.def("pixel_grid", &Dispersion::pixel_grid)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<Dispersion> >::*pbt1)(
        const std::vector<boost::shared_ptr<Dispersion> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<Dispersion> >, VectorDispersion)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<Dispersion> >::push_back))
REGISTER_LUA_END()

#endif
