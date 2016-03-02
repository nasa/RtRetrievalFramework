// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "lua_callback.h"
%}
%import "lua_state.i"
%fp_shared_ptr(FullPhysics::LuaCallback);
namespace FullPhysics {
%feature("director") LuaCallback;
class LuabindObject;
class LuaCallback {
public:
  LuaCallback(const LuaState& Ls);
  virtual ~LuaCallback();
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
       const boost::shared_ptr<LuabindObject>& Obj10
       ) = 0;
  std::string print_to_string() const;
};

}

%pythoncode %{
class LuaCallbackWrap(LuaCallback):
    def __init__(self,ls,f):
       LuaCallback.__init__(self, ls)
       self.ls = ls
       self.f = f
       t = inspect.getargspec(f)[0]
       self.narg = len(t)
       if(self.narg > 0 and t[0] == 'self'):
         self.narg = self.narg - 1
       if(self.narg > 10):
         raise "We only support up to 10 argument in a function. This limitation comes from luabind"

    def call(self,obj1, obj2, obj3, obj4, obj5, obj6, obj7, obj8, obj9,
             obj10):
       if(self.narg == 0):
         t = self.f()
       elif(self.narg == 1):
         t = self.f(obj1.cast_type())
       elif(self.narg == 2):
         t = self.f(obj1.cast_type(), obj2.cast_type())
       elif(self.narg == 3):
         t = self.f(obj1.cast_type(), obj2.cast_type(), obj3.cast_type())
       elif(self.narg == 4):
         t = self.f(obj1.cast_type(), obj2.cast_type(), obj3.cast_type(),
                    obj4.cast_type())
       elif(self.narg == 5):
         t = self.f(obj1.cast_type(), obj2.cast_type(), obj3.cast_type(),
                    obj4.cast_type(), obj5.cast_type())
       elif(self.narg == 6):
         t = self.f(obj1.cast_type(), obj2.cast_type(), obj3.cast_type(),
                    obj4.cast_type(), obj5.cast_type(), oj6.cast_type())
       elif(self.narg == 7):
         t = self.f(obj1.cast_type(), obj2.cast_type(), obj3.cast_type(),
                    obj4.cast_type(), obj5.cast_type(), obj6.cast_type(),
                    obj7.cast_type())
       elif(self.narg == 8):
         t = self.f(obj1.cast_type(), obj2.cast_type(), obj3.cast_type(),
                    obj4.cast_type(), obj5.cast_type(), obj6.cast_type(),
                    obj7.cast_type(), obj8.cast_type())
       elif(self.narg == 9):
         t = self.f(obj1.cast_type(), obj2.cast_type(), obj3.cast_type(),
                    obj4.cast_type(), obj5.cast_type(), obj6.cast_type(),
                    obj7.cast_type(), obj8.cast_type(), obj9.cast_type())
       elif(self.narg == 10):
         t = self.f(obj1.cast_type(), obj2.cast_type(), obj3.cast_type(),
                    obj4.cast_type(), obj5.cast_type(), obj6.cast_type(),
                    obj7.cast_type(), obj8.cast_type(), obj9.cast_type().
                    obj10.cast_type())
       else:
         raise "This shouldn't be able to happen"

       if(t == None):
           return LuabindObject.nil(self.ls)

       elif(isinstance(t, dict)):
           new_table = LuabindObject.nil(self.ls).new_table()
           for t_key, t_val in t.items():
              new_table[t_key] = t_val
           return LuabindObject(self.ls, new_table)

       elif(isinstance(t, list) or isinstance(t, tuple)):
           new_table = LuabindObject.nil(self.ls).new_table()
           for t_idx, t_val in enumerate(t):
               new_table[t_idx+1] = t_val
           return LuabindObject(self.ls, new_table)

       else:
         return LuabindObject(self.ls, t)
%}
