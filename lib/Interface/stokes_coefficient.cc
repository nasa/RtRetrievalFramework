#include "stokes_coefficient.h"

using namespace FullPhysics;
#ifdef HAVE_LUA
#include "register_lua.h"
blitz::Array<double, 2> stokes_coefficient_value(const StokesCoefficient& S)
{
  return S.stokes_coefficient().value();
}
REGISTER_LUA_CLASS(StokesCoefficient)
.def("stokes_coefficient", &StokesCoefficient::stokes_coefficient)
.def("stokes_coefficient_value", &stokes_coefficient_value)
REGISTER_LUA_END()
#endif

