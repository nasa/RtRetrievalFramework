#ifndef LUA_STATE_H
#define LUA_STATE_H
#include "printable.h"
#include "register_lua.h"
#include <boost/noncopyable.hpp>

namespace FullPhysics {
class LuabindObject;
// Actual class that holds the implementation. We use this level of
// indirection so we can make sure to maintain the lifetime of the
// underlying lua_State apart from the LuaState class.
class LuaStateImp : public boost::noncopyable {
public:
  LuaStateImp(lua_State* Ls, const std::string& Dir_name)
    : ls(Ls), dir_name(Dir_name) {}
  lua_State *ls;
  std::string dir_name;
  ~LuaStateImp();
};

/****************************************************************//**
  This is a light wrapper around the lua_state object. This maintains
  the lifetime of this object, as well as other house keeping chores.
*******************************************************************/
class LuaState : public Printable<LuaState> {
public:
  LuaState(const std::string& Dir_name = "./");
  virtual ~LuaState() {}
  void print(std::ostream& Os) const { Os << "LuaState";}

  static boost::shared_ptr<LuaState> load_file(const std::string& Fname);
  void do_file(const std::string& Fname);
  void run(const std::string& S);
  LuabindObject globals();
  LuabindObject registry();

//-----------------------------------------------------------------------
/// Base directory that we use when running Lua.
//-----------------------------------------------------------------------

  const std::string& base_dir_name() const {return lsimp->dir_name; }

//-----------------------------------------------------------------------
/// Lua state pointer. You can use this to do things with luabind (or
/// directly with Lua) that this class doesn't already support.
//-----------------------------------------------------------------------
  lua_State* lua_state() {return lsimp->ls;}
private:
  boost::shared_ptr<LuaStateImp> lsimp;
  LuaState(const boost::shared_ptr<LuaStateImp>& Lsimp)
    : lsimp(Lsimp) {}
};
}

// This needs to get defined after LuaState.
#include "luabind_object.h"
#endif
