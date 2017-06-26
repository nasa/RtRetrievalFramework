#include "spectrum.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
typedef const SpectralDomain& (Spectrum::*f1)(void) const;
typedef const SpectralRange& (Spectrum::*f2)(void) const;
REGISTER_LUA_CLASS(Spectrum)
.def("spectral_domain", ((f1) &Spectrum::spectral_domain))
.def("spectral_range", ((f2) &Spectrum::spectral_range))
REGISTER_LUA_END()
#endif
