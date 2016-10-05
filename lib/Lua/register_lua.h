#ifndef REGISTER_LUA_H
#define REGISTER_LUA_H

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}
#include <luabind/luabind.hpp>
#include <luabind/tag_function.hpp>
#include <boost/function.hpp>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "lua_state.h"

namespace FullPhysics {
/****************************************************************//**
  This class handles the registration of luabind class wrappers with
  Lua.

  One approach to this is to have a central function that registers
  everything, and as we add classes update that central function. An
  alternative is the one selected here, were we have a more
  decentralized registration. Classes set up the registration in their
  own area, and then simple get listed and needing registration in the
  file "register_lua.cc".  It would be nice to decentralize this
  completely, but I could never figure out a way to actually do this.

  So registration involves 2 steps:

  1. Add the registration code to the class code (e.g., for Foo, this
     is the file "foo.cc").
  2. Add the class to the list of classes in the function
     RegisterLua::register_lua found at "lib/Lua/register_lua.cc" 

  The registration code is cookie cutter, so we have macros to help
  do this. The registration is different depending on if we have a
  derived class with a base class, or a class that doesn't derive from
  another (or at least one that we want to tell Lua about).

  An example of this:

  In level_1b.cc:
  \code
   #ifdef HAVE_LUA
   #include "register_lua.h"
   REGISTER_LUA_CLASS(Level1b)
   REGISTER_LUA_END()
   #endif
  \endcode

  In level_1b_hdf.cc:
  \code
   #ifdef HAVE_LUA
   #include "register_lua.h"
   REGISTER_LUA_DERIVED_CLASS(Level1bAcos, Level1b)
     .def(luabind::constructor<std::string, std::string>())
   REGISTER_LUA_END()
   #endif
   \endcode

   Then in register_lua.cc, we add
   \code
   REGISTER_LUA_LIST(Level1b);
   REGISTER_LUA_LIST(Level1Hdf);
   \endcode

   Note that you don't need to put all the member functions into the
   Lua registration, just the ones you want to call in Lua. For many
   classes, this will just be the constructors. We use Lua
   configuration files for creating the objects needed in Level 2 Full
   physics, not to do major computation with it. That is more what we
   do with the Python wrappers. Lua is a small language that is ideal
   for integration in the C++ code, but it is no replacement for
   Python, nor is it meant to be.

   Pretty much all our classes are Printable. We've put the magic 
   incantation in place for classes in the macros (this ties the Lua
   function __tostring to the C++ code print_to_string). If you have 
   a class that is not printable, we'll need to add a macro to support that.

   We normally use Lua through our C++ code. It can be useful, 
   particularly when testing, to go the other way. We define the
   function "luaopen_fullphysics" to go the other way, call in 
   Lua like:

   \code
   require("libfull_physics")
   \endcode
   
   Note that you should use the installed library, like we do with 
   python (i.e., do a "make install").

   You will need to make sure that the library is on the PATH. Lua
   uses an odd syntax for its path, an example of using it would be

   \code
   LUA_CPATH=install/lib/?.so lua
   require "libfull_physics"
   l1b = Level1bAcos("filename","soundingid")
   \endcode
*******************************************************************/

class RegisterLua {
public:
  static int add_file_and_line(lua_State *ls);
  static void register_lua(lua_State *ls);
};
}

#define REGISTER_LUA_LIST(X) \
void register_lua_##X(lua_State*); register_lua_##X(ls)


#define REGISTER_LUA_CLASS(X) \
namespace FullPhysics { \
void register_lua_##X(lua_State *ls) { \
 LuabindObject::base_bridge.push_back(boost::shared_ptr<BridgeLuabindAndGenericBase>(new BridgeLuabindAndGeneric< X >)); \
luabind::module(ls) [ luabind::class_< X,boost::shared_ptr< X > >(#X) \
.def("__tostring", &X::print_to_string)

#define REGISTER_LUA_CLASS_NAME_WITH_GENERIC_OBJECT_BASE(X, Y) \
namespace FullPhysics { \
void register_lua_##Y(lua_State *ls) { \
 LuabindObject::base_bridge.push_back(boost::shared_ptr<BridgeLuabindAndGenericBase>(new BridgeLuabindAndGeneric< X >)); \
luabind::module(ls) [ luabind::class_< X,boost::shared_ptr< X > >(#Y)

#define REGISTER_LUA_CLASS_NAME(X, Y)		\
namespace FullPhysics { \
void register_lua_##Y(lua_State *ls) { \
luabind::module(ls) [ luabind::class_< X,boost::shared_ptr< X > >(#Y)

#define REGISTER_LUA_END() ]; } }

#define REGISTER_LUA_DERIVED_CLASS(X, Y) \
namespace FullPhysics { \
void register_lua_##X(lua_State *ls) { \
 LuabindObject::base_bridge.push_back(boost::shared_ptr<BridgeLuabindAndGenericBase>(new BridgeLuabindAndGeneric< Y >)); \
luabind::module(ls) [ luabind::class_<X,Y,boost::shared_ptr<Y> >(#X)

#endif
