#include "instrument_correction.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(InstrumentCorrection)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<InstrumentCorrection> >::*pbt1)(
        const std::vector<boost::shared_ptr<InstrumentCorrection> >::value_type&);
REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<InstrumentCorrection> >, VectorInstrumentCorrection)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<InstrumentCorrection> >::push_back))
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >::*pbt2)(
        const std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >::value_type&);
REGISTER_LUA_CLASS_NAME(std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >, VectorVectorInstrumentCorrection)
.def(luabind::constructor<>())
.def("push_back", ((pbt2) &std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >::push_back))
REGISTER_LUA_END()
#endif
