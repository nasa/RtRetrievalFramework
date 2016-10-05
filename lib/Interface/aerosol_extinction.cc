#include "aerosol_extinction.h"
#include <vector>
#include <boost/shared_ptr.hpp>

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(AerosolExtinction)
REGISTER_LUA_END()

typedef std::vector<boost::shared_ptr<AerosolExtinction> >::reference 
(std::vector<boost::shared_ptr<AerosolExtinction> >::*vaevt)(std::vector<boost::shared_ptr<AerosolExtinction> >::size_type);
REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<AerosolExtinction> >,
			VectorAerosolExtinction)
.def(luabind::constructor<>())
.def("size", &std::vector<boost::shared_ptr<AerosolExtinction> >::size)
.def("push_back", &std::vector<boost::shared_ptr<AerosolExtinction> >::push_back)
.def("value", ((vaevt) &std::vector<boost::shared_ptr<AerosolExtinction> >::operator[]))
REGISTER_LUA_END()
#endif
