#include "absorber.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Absorber)
.def("gas_index", &Absorber::gas_index)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Timer for optical_depth_each_layer.
//-----------------------------------------------------------------------

AccumulatedTimer Absorber::timer("Absorber optical_depth_each_layer");

//-----------------------------------------------------------------------
/// Map a gas name to the index number it appears in
/// optical_depth_each_layer. This return -1 if the Name is not one of
/// the gases.
//-----------------------------------------------------------------------

int Absorber::gas_index(const std::string& Name) const
{
  for(int i = 0; i < number_species(); ++i)
    if(gas_name(i) == Name)
      return i;
  return -1;
}


