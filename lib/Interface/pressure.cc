#include "pressure.h"

using namespace FullPhysics;
#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Pressure)
.def("max_number_level", &Pressure::max_number_level)
.def("pressure_grid", &Pressure::pressure_grid)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Return surface pressure, which is just the pressure at the bottom
/// level of pressure_grid.
///
/// This is in Pascals.
//-----------------------------------------------------------------------

AutoDerivativeWithUnit<double> Pressure::surface_pressure() const
{
  ArrayAdWithUnit<double, 1> p(pressure_grid());
  return p(p.rows() - 1);
}
