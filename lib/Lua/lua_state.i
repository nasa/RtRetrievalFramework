// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "lua_state.h"
%}

%fp_shared_ptr(FullPhysics::LuaState);
%fp_shared_ptr(FullPhysics::LuabindObject);

namespace FullPhysics {
class LuabindObject;
class LuaState {
public:
  LuaState(const std::string& Dir_name = "./");
  std::string print_to_string() const;

  static boost::shared_ptr<LuaState> load_file(const std::string& Fname);
  void do_file(const std::string& Fname);
  void run(const std::string& S);
  %python_attribute_nonconst(globals, LuabindObject)
  %python_attribute_nonconst(registry, LuabindObject)
};
}

