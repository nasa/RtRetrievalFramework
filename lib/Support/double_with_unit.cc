#include "double_with_unit.h"
#include "spectral_domain.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"

std::string double_with_unit_unit_get(const DoubleWithUnit& V)
{
  return V.units.name();
}

void double_with_unit_unit_set(DoubleWithUnit& V, std::string& Unit_name)
{
  V.units = Unit(Unit_name);
}

REGISTER_LUA_CLASS(DoubleWithUnit)
.def(luabind::constructor<double, const std::string&>())
.def_readwrite("value", &DoubleWithUnit::value)
.property("units", 
          &double_with_unit_unit_get,
          &double_with_unit_unit_set)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Variation of convert_wave that also handles the units of
/// sample_index. 
//-----------------------------------------------------------------------
    
DoubleWithUnit DoubleWithUnit::convert_wave
(const Unit& R, 
 const SpectralDomain& Pixel_grid) const
{
  if(units.is_commensurate(units::sample_index)) {
    int ind = (int) round(value) - 1;
    range_check(ind, 0, Pixel_grid.data().rows());
    DoubleWithUnit d(Pixel_grid.data()(ind), Pixel_grid.units());
    return d.convert_wave(R);
  } else {
    return convert_wave(R);
  }
}
