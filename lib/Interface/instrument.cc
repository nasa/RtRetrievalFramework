#include "instrument.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Instrument)
.def("number_spectrometer", &Instrument::number_spectrometer)
REGISTER_LUA_END()
#endif
