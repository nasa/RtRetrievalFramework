#include "lua_state.h"
#include "luabind_object.h"
#include "fp_exception.h"
#include "dir_change.h"
#include "fe_disable_exception.h"
using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Globals table.
//-----------------------------------------------------------------------

LuabindObject LuaState::globals() 
{
  luabind::object obj = luabind::globals(lua_state());
  boost::shared_ptr<LuaState> ls(new LuaState(lsimp));
  return LuabindObject(obj, ls); 
}

//-----------------------------------------------------------------------
/// Registery table.
//-----------------------------------------------------------------------

LuabindObject LuaState::registry()
{
  luabind::object obj = luabind::registry(lua_state());
  boost::shared_ptr<LuaState> ls(new LuaState(lsimp));
  return LuabindObject(obj, ls); 
}

//-----------------------------------------------------------------------
/// Create a new LuaState, and then open the given file. We run this
/// in the same directory as the given file, so for example a
/// configuration file can use relative paths for various files that
/// will be relative to the location of the configuration file, not
/// the directory we happen to be running in.
//-----------------------------------------------------------------------

boost::shared_ptr<LuaState> 
LuaState::load_file(const std::string& Fname)
{
  // Lua may cause floating point exceptions when loading the
  // configuration file. This is because it may copy garbage value,
  // which are never used. By chance, the garbage values may cause a
  // overflow. We suspend floating point exceptions when loading. 
  FeDisableException disable_fp;
  size_t t = Fname.find_last_of("/");
  std::string dirbase;
  std::string bname;
  if(t != std::string::npos) {
    dirbase = Fname.substr(0, t);
    bname = Fname.substr(t + 1);
  } else {
    dirbase = ".";
    bname = Fname;
  }
  dirbase += "/";
  boost::shared_ptr<LuaState> ls(new LuaState(dirbase));
  ls->do_file(bname);
  return ls;
}

//-----------------------------------------------------------------------
/// Create a new instance of Lua. Because it is often convenient to do
/// so, you can optionally give a directory that the we will change to
/// before running Lua code. This allows things like configuration
/// file to use relative paths from where the configuration file is
/// located rather than where we are running from.
//-----------------------------------------------------------------------

LuaState::LuaState(const std::string& Dir_name)
: lsimp(new LuaStateImp(luaL_newstate(), Dir_name))
{
  
  DirChange d(base_dir_name());
  luaL_openlibs(lua_state()); 
  luabind::open(lua_state());  
  RegisterLua::register_lua(lua_state());
}

LuaStateImp::~LuaStateImp()
{
  lua_close(ls);
  ls = 0;
}

//-----------------------------------------------------------------------
/// Load the given file and execute it.
//-----------------------------------------------------------------------

void LuaState::do_file(const std::string& Fname)
{
  DirChange d(base_dir_name());
  int status = luaL_dofile(lua_state(), Fname.c_str());
  if(status != 0) {
    // If we are here, then there is an error message on the stack
    Exception e;
    e << "Lua error: " << lua_tostring(lua_state(), -1) << "\n";
    lua_pop(lua_state(), 1);		// Remove error message from stack.
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Run the given Lua code.
//-----------------------------------------------------------------------

void LuaState::run(const std::string& S)
{
  DirChange d(base_dir_name());
  int status = luaL_dostring(lua_state(), S.c_str());
  if(status != 0) {
    // If we are here, then there is an error message on the stack
    Exception e;
    e << "Lua error: " << lua_tostring(lua_state(), -1) << "\n";
    lua_pop(lua_state(), 1);		// Remove error message from stack.
    throw e;
  }
}
