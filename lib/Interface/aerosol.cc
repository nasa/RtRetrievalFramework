#include "aerosol.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Aerosol)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Timer for Aerosol
//-----------------------------------------------------------------------

AccumulatedTimer Aerosol::timer("Aerosol");

