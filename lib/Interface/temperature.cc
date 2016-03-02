#include "temperature.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Temperature)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Return temperature at the pressure grid.
//-----------------------------------------------------------------------

ArrayAdWithUnit<double, 1> Temperature::temperature_grid(const Pressure& P) 
const
{
  ArrayAdWithUnit<double, 1> pgrid = P.pressure_grid();
  blitz::Array<AutoDerivative<double>, 1> res(pgrid.rows());
  Unit u = temperature(pgrid(0)).units;
  for(int i = 0; i < res.rows(); ++i)
    res(i) = temperature(pgrid(i)).convert(u).value;
  return ArrayAdWithUnit<double, 1>(ArrayAd<double, 1>(res), u);
}
