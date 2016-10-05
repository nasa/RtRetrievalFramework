#include "old_constant.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"

using namespace luabind;
namespace FullPhysics { 
  void register_lua_OldConstant(lua_State *ls) {
    // Try and retrieve global table Constants,
    // if not create it
    object const_table = globals(ls)["Constants"];

    if (type(const_table) == LUA_TNIL) {
      const_table = newtable(ls);
      globals(ls)["Constants"] = const_table;
    }

    // Add our constants to the global table
    const_table["speed_of_light"] = OldConstant::speed_of_light.value;
  } 
}
#endif
