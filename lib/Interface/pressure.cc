#include "pressure.h"

using namespace FullPhysics;
#ifdef HAVE_LUA
#include "register_lua.h"


blitz::Array<double, 1> pressure_level(const boost::shared_ptr<Pressure>& press) {
    return press->pressure_grid().value.value();
}

REGISTER_LUA_CLASS(Pressure)
.def("max_number_level", &Pressure::max_number_level)
.def("pressure_level", &pressure_level)
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
