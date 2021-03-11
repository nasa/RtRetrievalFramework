#include "ils_function.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_CLASS(IlsFunction)
.def("hdf_band_name", &IlsFunction::hdf_band_name)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<IlsFunction> >::*pbt1)(
        const std::vector<boost::shared_ptr<IlsFunction> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<IlsFunction> >, VectorIlsFunction)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<IlsFunction> >::push_back))
REGISTER_LUA_END()

#endif
