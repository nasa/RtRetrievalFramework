#include "ils_function.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(IlsFunction)
.def("hdf_band_name", &IlsFunction::hdf_band_name)
REGISTER_LUA_END()
#endif
