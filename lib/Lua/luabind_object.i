// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "luabind_object.h"
#include "altitude.h"
#include "state_vector.h"
#include "pressure.h"
#include "temperature.h"
%}
%import "generic_object.i"
%import "altitude.i"
%import "lua_state.i"
 // Really do want %include here. We have these pieces of code tightly
 // coupled, and easier if we include this in one wrapper
%include "lua_callback.i"
// This appears in LuaState, which is tightly coupled to this class.
// %shared_ptr(FullPhysics::LuabindObject);

// We get false warning for various version Array shadowing each other. This
// isn't a real problem, it has to do with how SWIG does the parsing (see
// http://sourceforge.net/p/swig/bugs/1053/). We just suppress this
%warnfilter(509) FullPhysics::LuabindObject;

%pythoncode %{
import inspect
%}

namespace FullPhysics {
class LuabindIterator {
public:
  LuabindIterator(const LuabindObject& Lbo,
		  const boost::shared_ptr<LuaState>& Ls);
  bool at_end() const;
  LuabindObject key() const;
  LuabindObject next();
};

class LuabindObject {
public:
  std::string print_to_string() const;
  bool is_nil() const;
  bool is_boolean() const;
  bool is_number() const;
  bool is_string() const;
  bool is_table() const;
  bool is_function() const;
  int length(int index = 0) const;
  LuabindObject call();
  LuabindObject call(const luabind::object& Arg1);
  LuabindObject call(const luabind::object& Arg1,
		     const luabind::object& Arg2);
  LuabindObject call(const luabind::object& Arg1,
		     const luabind::object& Arg2,
		     const luabind::object& Arg3);
  LuabindObject call(const luabind::object& Arg1,
		     const luabind::object& Arg2,
		     const luabind::object& Arg3,
		     const luabind::object& Arg4);
  LuabindObject call(const luabind::object& Arg1,
		     const luabind::object& Arg2,
		     const luabind::object& Arg3,
		     const luabind::object& Arg4,
		     const luabind::object& Arg5);
  LuabindObject call(const luabind::object& Arg1,
		     const luabind::object& Arg2,
		     const luabind::object& Arg3,
		     const luabind::object& Arg4,
		     const luabind::object& Arg5,
		     const luabind::object& Arg6);
  LuabindObject call(const luabind::object& Arg1,
		     const luabind::object& Arg2,
		     const luabind::object& Arg3,
		     const luabind::object& Arg4,
		     const luabind::object& Arg5,
		     const luabind::object& Arg6,
		     const luabind::object& Arg7);
  LuabindObject call(const luabind::object& Arg1,
		     const luabind::object& Arg2,
		     const luabind::object& Arg3,
		     const luabind::object& Arg4,
		     const luabind::object& Arg5,
		     const luabind::object& Arg6,
		     const luabind::object& Arg7,
		     const luabind::object& Arg8);
  const boost::shared_ptr<LuaState>& lua_state() const;
  luabind::object& object();
  LuabindObject new_table();
  static LuabindObject nil(const boost::shared_ptr<LuaState>& Ls);
  boost::shared_ptr<GenericObject> value_generic_object() const;
  void set_value(const std::string& Vname, 
		 const boost::shared_ptr<GenericObject>& V);
  void set_index(int Vidx, 
		 const boost::shared_ptr<GenericObject>& V);
  %extend {
    bool isvector_altitude () const
    { return $self->is_type_ptr<std::vector<boost::shared_ptr<FullPhysics::Altitude> > >(); }
    std::vector<boost::shared_ptr<FullPhysics::Altitude> >
      getvector_altitude() const
    { return *$self->value_ptr<std::vector<boost::shared_ptr<FullPhysics::Altitude> > >(); }

    bool isblitzarray1d_double () const
    { return $self->is_type<blitz::Array<double, 1> >(); }
    blitz::Array<double, 1> getblitzarray1d_double() const 
    { return $self->value<blitz::Array<double, 1> >(); }

    bool isblitzarray2d_double () const
    { return $self->is_type<blitz::Array<double, 2> >(); }
    blitz::Array<double, 2> getblitzarray2d_double() const 
    { return $self->value<blitz::Array<double, 2> >(); }

    // For passing retrieval flags
    bool isblitzarray1d_bool () const
    { return $self->is_type<blitz::Array<bool, 1> >(); }
    blitz::Array<bool, 1> getblitzarray1d_bool() const 
    { return $self->value<blitz::Array<bool, 1> >(); }

  }
 void set_value(const std::string& Vname,
	        const boost::shared_ptr<FullPhysics::LuaCallback>& V);
 void set_index(int Vidx,
	        const boost::shared_ptr<FullPhysics::LuaCallback>& V);

 void set_value(const std::string& Vname, const std::string& V);
 void set_value(const std::string& Vname, int V);
 void set_value(const std::string& Vname, double V);
 void set_value(const std::string& Vname, const LuabindObject& V);
 void set_value(const std::string& Vname, const blitz::Array<double, 1>& V);
 void set_value(const std::string& Vname, const blitz::Array<double, 2>& V);
 void set_value(const std::string& Vname, const blitz::Array<double, 3>& V);
 void set_value(const std::string& Vname, const blitz::Array<bool, 1>& V);
 void set_value(const std::string& Vname, const blitz::Array<bool, 2>& V);
 void set_value(const std::string& Vname, const blitz::Array<bool, 3>& V);
 void set_value(const std::string& Vname, const blitz::Array<int, 1>& V);
 void set_value(const std::string& Vname, const blitz::Array<int, 2>& V);
 void set_value(const std::string& Vname, const blitz::Array<int, 3>& V);

 void set_index(int Vidx, const std::string& V);
 void set_index(int Vidx, int V);
 void set_index(int Vidx, double V);
 void set_index(int Vidx, const LuabindObject& V);

 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const std::string& V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, int V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, double V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const LuabindObject& V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const blitz::Array<double, 1>& V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const blitz::Array<double, 2>& V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const blitz::Array<double, 3>& V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const blitz::Array<bool, 1>& V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const blitz::Array<bool, 2>& V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const blitz::Array<bool, 3>& V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const blitz::Array<int, 1>& V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const blitz::Array<int, 2>& V);
 LuabindObject(const boost::shared_ptr<LuaState>& Ls, const blitz::Array<int, 3>& V);
%extend {
  LuabindObject(const boost::shared_ptr<LuaState>& Ls, 
		const boost::shared_ptr<GenericObject>& Obj)
    {
      boost::shared_ptr<FullPhysics::LuabindObject> lb = FullPhysics::LuabindObject::create_luabind_object(Ls, Obj);
      return new FullPhysics::LuabindObject(*lb);
    }
  void set_value_bool(const std::string& Vname, bool V)
  { $self->set_value(Vname, V); }
  void set_index_bool(int Vidx, bool V)
  { $self->set_index(Vidx, V); }

  LuabindObject get_string(const std::string& Vname)
  { return (*$self)[Vname]; }
  LuabindObject get_index(int Vidx)
  { return $self->get_index(Vidx); }

  bool getbool() { return $self->value<bool>(); }
  double getnum() { return $self->value<double>(); }
  std::string getstring() { return $self->value<std::string>(); }
}
%pythoncode %{
def cast_type(self):
    '''Cast to underlying type'''
    t = "Types"
    if(self.is_nil()):
        return None
    v = self.value_generic_object()
    if(v is not None):
        return v
    elif(self.is_boolean()):
        return self.getbool()
    elif(self.is_number()):
        return self.getnum()
    elif(self.is_string()):
        return self.getstring()
    elif(self.is_table()):
        return self
    elif(self.is_function()):
        return self
    elif(self.isblitzarray1d_double()):
        return self.getblitzarray1d_double()
    elif(self.isblitzarray2d_double()):
        return self.getblitzarray2d_double()
    elif(self.isblitzarray1d_bool()):
        return self.getblitzarray1d_bool()
    elif(self.isvector_altitude()):
        return self.getvector_altitude()
    return self
    
def __getitem__(self, key):
    if isinstance(key, int):
        return self.get_index(key).cast_type()
    else:
        return self.get_string(key).cast_type()

def __setitem__(self, key, v):
    # Set up which set method to use
    if(isinstance(key, int)):
        if key <= 0:
            raise ValueError("Lua integer indexing is 1 based, not 0 based")
        set_method = self.set_index
    else:
        set_method = self.set_value

    # Perform any value specific work
    if(not isinstance(v, LuabindObject) and hasattr(v, '__call__')):
        set_method(key, LuaCallbackWrap(self.lua_state(), v))

    elif(isinstance(v, bool)):
        if(isinstance(key, int)):
            set_method = self.set_index_bool
        else:
            set_method = self.set_value_bool
        set_method(key, v)

    elif(isinstance(v, dict)):
        new_table = self.new_table()
        for v_key, v_val in v.items():
           new_table[v_key] = v_val
        set_method(key, new_table)

    elif(isinstance(v, list) or isinstance(v, tuple)):
        # Only support assigining lists or tuples for now,
        # since so many other objects are iterable but
        # need to be set differently, for instance LuabindObject
        # itself is iterable
        new_table = self.new_table()
        for v_idx, v_val in enumerate(v):
           new_table[v_idx+1] = v_val
        set_method(key, new_table)
    else:
        set_method(key, v)

def __getattr__(self, key):
    return self[key]

def __setattr__(self, key, v):
    if(key == "this"):
        self.__dict__[key] = v
    else:
        self[key] = v

%}
%pythoncode %{
def __call__(self, *args):
    luargs = []
    luargs2 = []
    for a in args:
        if(isinstance(a, LuabindObject)):
           t = a
        else:
           t = LuabindObject(self.lua_state(), a)
        luargs.append(t)
        luargs2.append(t.object())
    return self.call(*luargs2).cast_type()

def __len__(self):
    return self.cast_type().length()

def __iter__(self):
    if(self.is_table()):
        def luabind_gen():
            lb_iter = LuabindIterator(self, self.lua_state())
            while not lb_iter.at_end():
                yield lb_iter.next().cast_type()
        return luabind_gen()
    else:
        raise TypeError("Can not iterate on non table LuabindObject")

def __dir__(self):
    if(self.is_table()):
        obj_contents = []
        lb_iter = LuabindIterator(self, self.lua_state())
        while not lb_iter.at_end():
            key_item = lb_iter.key().cast_type()
            # Only add strings to dir() listing,
            # because in Lua arrays would have numeric
            # keys and these are not accessible through
            # getattr
            if isinstance(key_item, str):
                obj_contents.append( key_item  )
            lb_iter.next()
        return obj_contents
    else:
        item = self.cast_type()
        # Try not to engage in infinite recursion!
        if isinstance(item, LuabindObject):
            return []
        else:
            return dir(item)
%}
};
}
