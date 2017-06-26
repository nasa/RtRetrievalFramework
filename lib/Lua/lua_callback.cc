#include "lua_callback.h"
using namespace FullPhysics;
#ifdef HAVE_LUA
#include "register_lua.h"
typedef luabind::object (LuaCallback::*f0)();
typedef luabind::object (LuaCallback::*f1)(const luabind::object&);
typedef luabind::object (LuaCallback::*f2)(const luabind::object&, 
					   const luabind::object&);
typedef luabind::object (LuaCallback::*f3)(const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&);
typedef luabind::object (LuaCallback::*f4)(const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&);
typedef luabind::object (LuaCallback::*f5)(const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&);
typedef luabind::object (LuaCallback::*f6)(const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&);
typedef luabind::object (LuaCallback::*f7)(const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&);
typedef luabind::object (LuaCallback::*f8)(const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&);
typedef luabind::object (LuaCallback::*f9)(const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&);
typedef luabind::object (LuaCallback::*f10)(const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&, 
					   const luabind::object&);
REGISTER_LUA_CLASS(LuaCallback)
.def("__call", (f0) &LuaCallback::__call)
.def("__call", (f1) &LuaCallback::__call)
.def("__call", (f2) &LuaCallback::__call)
.def("__call", (f3) &LuaCallback::__call)
.def("__call", (f4) &LuaCallback::__call)
.def("__call", (f5) &LuaCallback::__call)
.def("__call", (f6) &LuaCallback::__call)
.def("__call", (f7) &LuaCallback::__call)
.def("__call", (f8) &LuaCallback::__call)
.def("__call", (f9) &LuaCallback::__call)
.def("__call", (f10) &LuaCallback::__call)
REGISTER_LUA_END()

#endif
