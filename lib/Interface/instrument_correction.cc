#include "instrument_correction.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(InstrumentCorrection)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<InstrumentCorrection> >,
			VectorInstrumentCorrection)
.def(luabind::constructor<>())
.def("push_back", &std::vector<boost::shared_ptr<InstrumentCorrection> >::push_back)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >,
			VectorVectorInstrumentCorrection)
.def(luabind::constructor<>())
.def("push_back", &std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >::push_back)
REGISTER_LUA_END()
#endif
