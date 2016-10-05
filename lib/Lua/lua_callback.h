#ifndef LUA_CALLBACK_H
#define LUA_CALLBACK_H
#include "printable.h"
#include "lua_state.h"

namespace FullPhysics {
/****************************************************************//**
  This is a simple object to call a callback that can be used in Lua.
  The function can then execute C++ or Python code.
*******************************************************************/
class LuaCallback : public Printable<LuaCallback> {
public:
  LuaCallback(const LuaState& Ls) : ls(new LuaState(Ls)) {}
  virtual ~LuaCallback() {}
  virtual void print(std::ostream& Os) const { Os << "LuaCallback"; }
  virtual boost::shared_ptr<LuabindObject> 
  call(const boost::shared_ptr<LuabindObject>& Obj1,
       const boost::shared_ptr<LuabindObject>& Obj2,
       const boost::shared_ptr<LuabindObject>& Obj3,
       const boost::shared_ptr<LuabindObject>& Obj4,
       const boost::shared_ptr<LuabindObject>& Obj5,
       const boost::shared_ptr<LuabindObject>& Obj6,
       const boost::shared_ptr<LuabindObject>& Obj7,
       const boost::shared_ptr<LuabindObject>& Obj8,
       const boost::shared_ptr<LuabindObject>& Obj9,
       const boost::shared_ptr<LuabindObject>& Obj10) = 0;
  luabind::object __call(const luabind::object& obj1, 
			 const luabind::object& obj2,
			 const luabind::object& obj3,
			 const luabind::object& obj4,
			 const luabind::object& obj5,
			 const luabind::object& obj6,
			 const luabind::object& obj7,
			 const luabind::object& obj8,
			 const luabind::object& obj9,
			 const luabind::object& obj10) 
  { boost::shared_ptr<LuabindObject> lobj1(new LuabindObject(obj1, ls)); 
    boost::shared_ptr<LuabindObject> lobj2(new LuabindObject(obj2, ls));
    boost::shared_ptr<LuabindObject> lobj3(new LuabindObject(obj3, ls));
    boost::shared_ptr<LuabindObject> lobj4(new LuabindObject(obj4, ls));
    boost::shared_ptr<LuabindObject> lobj5(new LuabindObject(obj5, ls));
    boost::shared_ptr<LuabindObject> lobj6(new LuabindObject(obj6, ls));
    boost::shared_ptr<LuabindObject> lobj7(new LuabindObject(obj7, ls));
    boost::shared_ptr<LuabindObject> lobj8(new LuabindObject(obj8, ls));
    boost::shared_ptr<LuabindObject> lobj9(new LuabindObject(obj9, ls));
    boost::shared_ptr<LuabindObject> lobj10(new LuabindObject(obj10, ls));
    return call(lobj1, lobj2, lobj3, lobj4, lobj5, lobj6,
		lobj7, lobj8, lobj9, lobj10)->object();
  }
  luabind::object __call(const luabind::object& obj1, 
			 const luabind::object& obj2,
			 const luabind::object& obj3,
			 const luabind::object& obj4,
			 const luabind::object& obj5,
			 const luabind::object& obj6,
			 const luabind::object& obj7,
			 const luabind::object& obj8,
			 const luabind::object& obj9)
  { boost::shared_ptr<LuabindObject> lobj1(new LuabindObject(obj1, ls)); 
    boost::shared_ptr<LuabindObject> lobj2(new LuabindObject(obj2, ls));
    boost::shared_ptr<LuabindObject> lobj3(new LuabindObject(obj3, ls));
    boost::shared_ptr<LuabindObject> lobj4(new LuabindObject(obj4, ls));
    boost::shared_ptr<LuabindObject> lobj5(new LuabindObject(obj5, ls));
    boost::shared_ptr<LuabindObject> lobj6(new LuabindObject(obj6, ls));
    boost::shared_ptr<LuabindObject> lobj7(new LuabindObject(obj7, ls));
    boost::shared_ptr<LuabindObject> lobj8(new LuabindObject(obj8, ls));
    boost::shared_ptr<LuabindObject> lobj9(new LuabindObject(obj9, ls));
    boost::shared_ptr<LuabindObject> nil(new LuabindObject(luabind::object(), ls));
    return call(lobj1, lobj2, lobj3, lobj4, lobj5, lobj6,
		lobj7, lobj8, lobj9, nil)->object();
  }
  luabind::object __call(const luabind::object& obj1, 
			 const luabind::object& obj2,
			 const luabind::object& obj3,
			 const luabind::object& obj4,
			 const luabind::object& obj5,
			 const luabind::object& obj6,
			 const luabind::object& obj7,
			 const luabind::object& obj8)
  { boost::shared_ptr<LuabindObject> lobj1(new LuabindObject(obj1, ls)); 
    boost::shared_ptr<LuabindObject> lobj2(new LuabindObject(obj2, ls));
    boost::shared_ptr<LuabindObject> lobj3(new LuabindObject(obj3, ls));
    boost::shared_ptr<LuabindObject> lobj4(new LuabindObject(obj4, ls));
    boost::shared_ptr<LuabindObject> lobj5(new LuabindObject(obj5, ls));
    boost::shared_ptr<LuabindObject> lobj6(new LuabindObject(obj6, ls));
    boost::shared_ptr<LuabindObject> lobj7(new LuabindObject(obj7, ls));
    boost::shared_ptr<LuabindObject> lobj8(new LuabindObject(obj8, ls));
    boost::shared_ptr<LuabindObject> nil(new LuabindObject(luabind::object(), ls));
    return call(lobj1, lobj2, lobj3, lobj4, lobj5, lobj6,
		lobj7, lobj8, nil, nil)->object();
  }
  luabind::object __call(const luabind::object& obj1, 
			 const luabind::object& obj2,
			 const luabind::object& obj3,
			 const luabind::object& obj4,
			 const luabind::object& obj5,
			 const luabind::object& obj6,
			 const luabind::object& obj7
			 )
  { boost::shared_ptr<LuabindObject> lobj1(new LuabindObject(obj1, ls)); 
    boost::shared_ptr<LuabindObject> lobj2(new LuabindObject(obj2, ls));
    boost::shared_ptr<LuabindObject> lobj3(new LuabindObject(obj3, ls));
    boost::shared_ptr<LuabindObject> lobj4(new LuabindObject(obj4, ls));
    boost::shared_ptr<LuabindObject> lobj5(new LuabindObject(obj5, ls));
    boost::shared_ptr<LuabindObject> lobj6(new LuabindObject(obj6, ls));
    boost::shared_ptr<LuabindObject> lobj7(new LuabindObject(obj7, ls));
    boost::shared_ptr<LuabindObject> nil(new LuabindObject(luabind::object(), ls));
    return call(lobj1, lobj2, lobj3, lobj4, lobj5, lobj6,
		lobj7, nil, nil, nil)->object();
  }
  luabind::object __call(const luabind::object& obj1, 
			 const luabind::object& obj2,
			 const luabind::object& obj3,
			 const luabind::object& obj4,
			 const luabind::object& obj5,
			 const luabind::object& obj6
			 )
  { boost::shared_ptr<LuabindObject> lobj1(new LuabindObject(obj1, ls)); 
    boost::shared_ptr<LuabindObject> lobj2(new LuabindObject(obj2, ls));
    boost::shared_ptr<LuabindObject> lobj3(new LuabindObject(obj3, ls));
    boost::shared_ptr<LuabindObject> lobj4(new LuabindObject(obj4, ls));
    boost::shared_ptr<LuabindObject> lobj5(new LuabindObject(obj5, ls));
    boost::shared_ptr<LuabindObject> lobj6(new LuabindObject(obj6, ls));
    boost::shared_ptr<LuabindObject> nil(new LuabindObject(luabind::object(), ls));
    return call(lobj1, lobj2, lobj3, lobj4, lobj5, lobj6,
		nil, nil, nil, nil)->object();
  }
  luabind::object __call(const luabind::object& obj1, 
			 const luabind::object& obj2,
			 const luabind::object& obj3,
			 const luabind::object& obj4,
			 const luabind::object& obj5
			 )
  { boost::shared_ptr<LuabindObject> lobj1(new LuabindObject(obj1, ls)); 
    boost::shared_ptr<LuabindObject> lobj2(new LuabindObject(obj2, ls));
    boost::shared_ptr<LuabindObject> lobj3(new LuabindObject(obj3, ls));
    boost::shared_ptr<LuabindObject> lobj4(new LuabindObject(obj4, ls));
    boost::shared_ptr<LuabindObject> lobj5(new LuabindObject(obj5, ls));
    boost::shared_ptr<LuabindObject> nil(new LuabindObject(luabind::object(), ls));
    return call(lobj1, lobj2, lobj3, lobj4, lobj5, nil,
		nil, nil, nil, nil)->object();
  }
  luabind::object __call(const luabind::object& obj1, 
			 const luabind::object& obj2,
			 const luabind::object& obj3,
			 const luabind::object& obj4
			 )
  { boost::shared_ptr<LuabindObject> lobj1(new LuabindObject(obj1, ls)); 
    boost::shared_ptr<LuabindObject> lobj2(new LuabindObject(obj2, ls));
    boost::shared_ptr<LuabindObject> lobj3(new LuabindObject(obj3, ls));
    boost::shared_ptr<LuabindObject> lobj4(new LuabindObject(obj4, ls));
    boost::shared_ptr<LuabindObject> nil(new LuabindObject(luabind::object(), ls));
    return call(lobj1, lobj2, lobj3, lobj4, nil, nil,
		nil, nil, nil, nil)->object();
  }
  luabind::object __call(const luabind::object& obj1, 
			 const luabind::object& obj2,
			 const luabind::object& obj3
			 )
  { boost::shared_ptr<LuabindObject> lobj1(new LuabindObject(obj1, ls)); 
    boost::shared_ptr<LuabindObject> lobj2(new LuabindObject(obj2, ls));
    boost::shared_ptr<LuabindObject> lobj3(new LuabindObject(obj3, ls));
    boost::shared_ptr<LuabindObject> nil(new LuabindObject(luabind::object(), ls));
    return call(lobj1, lobj2, lobj3, nil, nil, nil,
		nil, nil, nil, nil)->object();
  }
  luabind::object __call(const luabind::object& obj1, 
			 const luabind::object& obj2
			 )
  { boost::shared_ptr<LuabindObject> lobj1(new LuabindObject(obj1, ls)); 
    boost::shared_ptr<LuabindObject> lobj2(new LuabindObject(obj2, ls));
    boost::shared_ptr<LuabindObject> nil(new LuabindObject(luabind::object(), ls));
    return call(lobj1, lobj2, nil, nil, nil, nil,
		nil, nil, nil, nil)->object();
  }
  luabind::object __call(const luabind::object& obj1
			 )
  { boost::shared_ptr<LuabindObject> lobj1(new LuabindObject(obj1, ls)); 
    boost::shared_ptr<LuabindObject> nil(new LuabindObject(luabind::object(), ls));
    return call(lobj1, nil, nil, nil, nil, nil,
		nil, nil, nil, nil)->object();
  }
  luabind::object __call()
  {
    boost::shared_ptr<LuabindObject> nil(new LuabindObject(luabind::object(), ls));
    return call(nil, nil, nil, nil, nil, nil,
		nil, nil, nil, nil)->object();
  }
private:
  boost::shared_ptr<LuaState> ls;
};

}
#endif

