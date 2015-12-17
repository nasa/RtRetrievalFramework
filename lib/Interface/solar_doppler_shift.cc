#include "solar_doppler_shift.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(SolarDopplerShift)
.def("solar_distance", &SolarDopplerShift::solar_distance)
REGISTER_LUA_END()
#endif
