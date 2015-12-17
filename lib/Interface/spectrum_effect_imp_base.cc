#include "spectrum_effect_imp_base.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(SpectrumEffectImpBase, SpectrumEffect)
REGISTER_LUA_END()
#endif
