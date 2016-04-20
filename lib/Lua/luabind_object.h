#ifndef LUABIND_OBJECT_H
#define LUABIND_OBJECT_H
#include "printable.h"
#include "register_lua.h"
#include "fp_exception.h"
#include "dir_change.h"
#include "lua_state.h"
#include <blitz/array.h>

namespace FullPhysics {
  class LuabindObject;

/****************************************************************//**
  For use with python, it is useful to map all the LuabindObject that
  are a class defined in full physics (or derived from those classes)
  to and from a GenericObject. This is pretty much pointless in C++,
  because you can't do anything with a GenericObject. But we have
  logic in the python wrappers to map a GenericObject to the most
  derived python object. This class is used to support this
  transformation in LuabindObject.
*******************************************************************/
class BridgeLuabindAndGenericBase {
public:
  virtual ~BridgeLuabindAndGenericBase() {}

//-----------------------------------------------------------------------
/// Create a LuabindObject from the given type, or return a null
/// pointer if we can't.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<LuabindObject> 
  luabind_object(const boost::shared_ptr<LuaState>& Ls,
		 const boost::shared_ptr<GenericObject>& Obj)
    const = 0;

//-----------------------------------------------------------------------
/// Return a GenericObject if we can convert, or a null pointer if
/// we can't.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<GenericObject> value(const LuabindObject& Lobj) 
    const = 0;

//-----------------------------------------------------------------------
/// Set a value if we can. Return true if we set the value, false otherwise
//-----------------------------------------------------------------------

  virtual bool set_value(LuabindObject& Lobj,
			 const std::string& Vname, 
			 const boost::shared_ptr<GenericObject>& Obj) = 0;

//-----------------------------------------------------------------------
/// Set a value if we can. Return true if we set the value, false otherwise
//-----------------------------------------------------------------------

  virtual bool set_index(LuabindObject& Lobj,
			 int Vidx, 
			 const boost::shared_ptr<GenericObject>& Obj) = 0;
};

/****************************************************************//**
  This is a light wrapper around the luabind::object. This adds a 
  few convenience routines.
*******************************************************************/
class LuabindObject : public Printable<LuabindObject> {
public:
  LuabindObject() {}

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
  LuabindObject(const luabind::object& Obj, 
		const boost::shared_ptr<LuaState>& Ls) : ls(Ls), obj(Obj) {}

//-----------------------------------------------------------------------
/// Return a nil value.
//-----------------------------------------------------------------------
  static LuabindObject nil(const boost::shared_ptr<LuaState>& Ls)
  { return LuabindObject(luabind::object(),Ls); }

//-----------------------------------------------------------------------
/// Conversion constructor
//-----------------------------------------------------------------------

  template<class T> LuabindObject(const boost::shared_ptr<LuaState>& Ls, 
				  const T& V) 
    : ls(Ls), obj(Ls->lua_state(), V) {}

// Nice not to have to treat LuabindObject as a special case, so just
// copy if T happens to be type LuabindObject
  LuabindObject(const boost::shared_ptr<LuaState>& Ls, 
		const LuabindObject& V) 
    : ls(V.lua_state()), obj(V.object()) {}

  virtual ~LuabindObject() {}

  void print(std::ostream& Os) const { Os << "LuabindObject";}

//-----------------------------------------------------------------------
/// Map used for mapping base classes to and from Lua and
/// GenericObject. See discussion above BridgeLuabindAndGenericBase.
//-----------------------------------------------------------------------

  static std::vector<boost::shared_ptr<BridgeLuabindAndGenericBase> > 
  base_bridge;

//-----------------------------------------------------------------------
/// Map used for mapping derived class to and from Lua and
/// GenericObject. See discussion above BridgeLuabindAndGenericBase.
//-----------------------------------------------------------------------

  static std::vector<boost::shared_ptr<BridgeLuabindAndGenericBase> > 
  derived_bridge;

  boost::shared_ptr<GenericObject> value_generic_object() const;
  void set_value(const std::string& Vname, 
				const boost::shared_ptr<GenericObject>& V);
  void set_index(int Vidx, const boost::shared_ptr<GenericObject>& V);
  static boost::shared_ptr<LuabindObject> 
  create_luabind_object(const boost::shared_ptr<LuaState>& Ls,
			const boost::shared_ptr<GenericObject>& Obj);

//-----------------------------------------------------------------------
/// Return a variable found in this table as another LuabindObject.
//-----------------------------------------------------------------------

  LuabindObject operator[](const std::string& Vname) const
  { luabind::object r = obj[Vname]; return LuabindObject(r, ls);}

  LuabindObject get_index(int Vidx) const
    { luabind::object r = obj[Vidx]; return LuabindObject(r, ls);}

//-----------------------------------------------------------------------
/// Return value of a variable.
//-----------------------------------------------------------------------

  template<class T> T
  value() const
  {
    try {
      return luabind::object_cast<T>(obj);
    } catch(std::exception eoriginal) {
      Exception e;
      e << "Error converting lua value to\n"
        << "Requested type: " << typeid(T).name() << "\n";
      if (is_nil())
        e << "Source object is nil";
      else
        e << "Source object is of type: " << typeid(obj).name();
      throw e;
    }
  }

//-----------------------------------------------------------------------
/// Return value of a variable found in this table. Shortcut of value
/// that lets you leave the boost::shared_ptr off in the type,
/// since this is so common. 
//-----------------------------------------------------------------------

  template<class T> boost::shared_ptr<T> 
  value_ptr() const
  {
    return value<boost::shared_ptr<T> >();
  }

//-----------------------------------------------------------------------
/// Test if value in a table is the given type. Return true if it is,
/// false otherwise.
//-----------------------------------------------------------------------

  template<class T> bool is_type() const
  { return luabind::object_cast_nothrow<T>(obj) != boost::none; }

//-----------------------------------------------------------------------
/// Test if value in a table is the given type. This is a shortcut for
/// is_type where you can leave off the boost::shared_ptr, since that
/// is so common.
//-----------------------------------------------------------------------

  template<class T> bool is_type_ptr() const
  {
    return is_type<boost::shared_ptr<T> >();
  }

//-----------------------------------------------------------------------
/// Set a value.
//-----------------------------------------------------------------------

  template<class T> void set_value(const std::string& Vname, 
				   const boost::shared_ptr<T>& V)
  {
    obj[Vname] = V;
  }
  template<class T, int D> void set_value(const std::string& Vname, 
					  const blitz::Array<T, D>& V)
  {
    obj[Vname] = V.copy();
  }
  void set_value(const std::string& Vname, const std::string& V)
  {
    obj[Vname] = V;
  }
  void set_value(const std::string& Vname, const char* V)
  {
    obj[Vname] = V;
  }
  void set_value(const std::string& Vname, int V)
  {
    obj[Vname] = V;
  }
  void set_value(const std::string& Vname, double V)
  {
    obj[Vname] = V;
  }
  void set_value(const std::string& Vname, bool V)
  {
    obj[Vname] = V;
  }
  void set_value(const std::string& Vname, const LuabindObject& V)
  {
    obj[Vname] = V.object();
  }

//-----------------------------------------------------------------------
/// Set an index
//-----------------------------------------------------------------------

  template<class T> void set_index(int Vidx, 
				   const boost::shared_ptr<T>& V)
  {
    obj[Vidx] = V;
  }
  void set_index(int Vidx, const std::string& V)
  {
    obj[Vidx] = V;
  }
  void set_index(int Vidx, const char* V)
  {
    obj[Vidx] = V;
  }
  void set_index(int Vidx, int V)
  {
    obj[Vidx] = V;
  }
  void set_index(int Vidx, double V)
  {
    obj[Vidx] = V;
  }
  void set_index(int Vidx, bool V)
  {
    obj[Vidx] = V;
  }
  void set_index(int Vidx, const LuabindObject& V)
  {
    obj[Vidx] = V.object();
  }

//-----------------------------------------------------------------------
/// Test type of an object.
//-----------------------------------------------------------------------

  bool is_nil() const { return luabind::type(obj) == LUA_TNIL; }
  bool is_boolean() const { return luabind::type(obj) == LUA_TBOOLEAN; }
  bool is_number() const { return luabind::type(obj) == LUA_TNUMBER; }
  bool is_string() const { return luabind::type(obj) == LUA_TSTRING; }
  bool is_table() const { return luabind::type(obj) == LUA_TTABLE; }
  bool is_function() const { return luabind::type(obj) == LUA_TFUNCTION; }

//-----------------------------------------------------------------------
/// Return the "length" of the value at the given index using the Lua
/// interpreter method that evaluates as the length operator
//-----------------------------------------------------------------------

  int length(int index = 0) const
  {
    return lua_rawlen(obj.interpreter(), index);
  }

//-----------------------------------------------------------------------
/// Call a function pointed to by this object.
//-----------------------------------------------------------------------
  LuabindObject call() 
  { 
    DirChange dc(ls->base_dir_name());
    try { 
      return LuabindObject((luabind::object) 
			   luabind::call_function<luabind::object>(obj), ls);
    } catch(const luabind::error& eoriginal) {
      Exception e;
      e << "Lua error: " << lua_tostring(lua_state()->lua_state(), -1) << "\n";
      lua_pop(lua_state()->lua_state(), 1); // Remove error message from stack.
      throw e;
    }
  }
  template<class T1>
  LuabindObject call(const T1& Arg1) 
  { 
    DirChange dc(ls->base_dir_name());
    try {
      return LuabindObject((luabind::object) 
			   luabind::call_function<luabind::object>(obj, Arg1), 
			   ls); 
    } catch(const luabind::error& eoriginal) {
      Exception e;
      e << "Lua error: " << lua_tostring(lua_state()->lua_state(), -1) << "\n";
      lua_pop(lua_state()->lua_state(), 1); // Remove error message from stack.
      throw e;
    }
  }
  template<class T1, class T2>
  LuabindObject call(const T1& Arg1, const T2& Arg2) 
  { 
    DirChange dc(ls->base_dir_name());
    try {
      return LuabindObject((luabind::object) 
			   luabind::call_function<luabind::object>
			   (obj, Arg1, Arg2), 
			   ls); 
    } catch(const luabind::error& eoriginal) {
      Exception e;
      e << "Lua error: " << lua_tostring(lua_state()->lua_state(), -1) << "\n";
      lua_pop(lua_state()->lua_state(), 1); // Remove error message from stack.
      throw e;
    }
  }
  template<class T1, class T2, class T3>
  LuabindObject call(const T1& Arg1, const T2& Arg2, const T3& Arg3) 
  { 
    DirChange dc(ls->base_dir_name());
    try {
      return LuabindObject((luabind::object) 
			   luabind::call_function<luabind::object>
			   (obj, Arg1, Arg2, Arg3), 
			   ls); 
    } catch(const luabind::error& eoriginal) {
      Exception e;
      e << "Lua error: " << lua_tostring(lua_state()->lua_state(), -1) << "\n";
      lua_pop(lua_state()->lua_state(), 1); // Remove error message from stack.
      throw e;
    }
  }
  template<class T1, class T2, class T3, class T4>
  LuabindObject call(const T1& Arg1, const T2& Arg2, const T3& Arg3,
		     const T4& Arg4) 
  { 
    DirChange dc(ls->base_dir_name());
    try {
      return LuabindObject((luabind::object) 
			   luabind::call_function<luabind::object>
			   (obj, Arg1, Arg2, Arg3, Arg4), 
			   ls); 
    } catch(const luabind::error& eoriginal) {
      Exception e;
      e << "Lua error: " << lua_tostring(lua_state()->lua_state(), -1) << "\n";
      lua_pop(lua_state()->lua_state(), 1); // Remove error message from stack.
      throw e;
    }

  }
  template<class T1, class T2, class T3, class T4, class T5>
  LuabindObject call(const T1& Arg1, const T2& Arg2, const T3& Arg3,
		     const T4& Arg4, const T5& Arg5) 
  { 
    DirChange dc(ls->base_dir_name());
    try{ 
      return LuabindObject((luabind::object) 
			   luabind::call_function<luabind::object>
			   (obj, Arg1, Arg2, Arg3, Arg4, Arg5), 
			   ls); 
    } catch(const luabind::error& eoriginal) {
      Exception e;
      e << "Lua error: " << lua_tostring(lua_state()->lua_state(), -1) << "\n";
      lua_pop(lua_state()->lua_state(), 1); // Remove error message from stack.
      throw e;
    }
  }
  template<class T1, class T2, class T3, class T4, class T5, class T6>
  LuabindObject call(const T1& Arg1, const T2& Arg2, const T3& Arg3,
		     const T4& Arg4, const T5& Arg5, const T6& Arg6) 
  { 
    DirChange dc(ls->base_dir_name());
    try { 
      return LuabindObject((luabind::object) 
			   luabind::call_function<luabind::object>
			   (obj, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6), 
			   ls); 
    } catch(const luabind::error& eoriginal) {
      Exception e;
      e << "Lua error: " << lua_tostring(lua_state()->lua_state(), -1) << "\n";
      lua_pop(lua_state()->lua_state(), 1); // Remove error message from stack.
      throw e;
    }
  }
  template<class T1, class T2, class T3, class T4, class T5, class T6,
	   class T7>
  LuabindObject call(const T1& Arg1, const T2& Arg2, const T3& Arg3,
		     const T4& Arg4, const T5& Arg5, const T6& Arg6,
		     const T7& Arg7) 
  { 
    DirChange dc(ls->base_dir_name());
    try {
      return LuabindObject((luabind::object) 
			   luabind::call_function<luabind::object>
			   (obj, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7), 
			   ls); 
    } catch(const luabind::error& eoriginal) {
      Exception e;
      e << "Lua error: " << lua_tostring(lua_state()->lua_state(), -1) << "\n";
      lua_pop(lua_state()->lua_state(), 1); // Remove error message from stack.
      throw e;
    }
  }
  template<class T1, class T2, class T3, class T4, class T5, class T6,
	   class T7, class T8>
  LuabindObject call(const T1& Arg1, const T2& Arg2, const T3& Arg3,
		     const T4& Arg4, const T5& Arg5, const T6& Arg6,
		     const T7& Arg7, const T8& Arg8) 
  { 
    DirChange dc(ls->base_dir_name());
    try {
      return LuabindObject((luabind::object) 
			   luabind::call_function<luabind::object>
			   (obj, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7, 
			    Arg8), ls); 
    } catch(const luabind::error& eoriginal) {
      Exception e;
      e << "Lua error: " << lua_tostring(lua_state()->lua_state(), -1) << "\n";
      lua_pop(lua_state()->lua_state(), 1); // Remove error message from stack.
      throw e;
    }
  }

//-----------------------------------------------------------------------
/// Creates a new table object
//-----------------------------------------------------------------------
  
  LuabindObject new_table() {
    return LuabindObject(luabind::newtable(ls->lua_state()), ls);
  }
  
//-----------------------------------------------------------------------
/// Underlying object. You can use this to do things with luabind that
/// this class doesn't already support.
//-----------------------------------------------------------------------

  luabind::object& object() {return obj;}
  const luabind::object& object() const {return obj;}

//-----------------------------------------------------------------------
/// Lua state pointer. You can use this to do things with luabind (or
/// directly with Lua) that this class doesn't already support.
//-----------------------------------------------------------------------

  const boost::shared_ptr<LuaState>& lua_state() const
  {return ls;}
private:
  // Note order here is important. We need obj to be deleted first,
  // before ls is.
  boost::shared_ptr<LuaState> ls;
  luabind::object obj;
};

/****************************************************************//**
  This is a light wrapper around the luabind::iterator. This adds a 
  few convenience routines to make usage in Python much easier.
*******************************************************************/

class LuabindIterator : public luabind::iterator {
public:
  //-----------------------------------------------------------------------
  /// Constructor.
  //-----------------------------------------------------------------------
  LuabindIterator(const LuabindObject& Lbo,
		  const boost::shared_ptr<LuaState>& Ls) : luabind::iterator(Lbo.object()), ls(Ls) {}

  bool at_end() const { return *this == luabind::iterator(); }
  LuabindObject key() const { return LuabindObject(luabind::iterator::key(), ls); }
  LuabindObject next() { LuabindObject res = LuabindObject(operator*(), ls); (*this)++; return res; }
private:
  boost::shared_ptr<LuaState> ls;
};

/****************************************************************//**
  Instance of BridgeLuabindAndGenericBase for a particular type
*******************************************************************/
template<class T>  class BridgeLuabindAndGeneric : 
    public BridgeLuabindAndGenericBase {
public:
  BridgeLuabindAndGeneric() {}
  virtual ~BridgeLuabindAndGeneric() {}
  virtual boost::shared_ptr<LuabindObject> 
  luabind_object(const boost::shared_ptr<LuaState>& Ls,
		 const boost::shared_ptr<GenericObject>& Obj)
    const 
  {
    boost::shared_ptr<T> v = boost::dynamic_pointer_cast<T>(Obj);
    boost::shared_ptr<LuabindObject> res;
    if(v)
      res.reset(new LuabindObject(Ls, v));
    return res;
  }
  
  virtual boost::shared_ptr<GenericObject> value(const LuabindObject& Lobj) 
    const
  {
    if(Lobj.is_type_ptr<T>())
      return boost::dynamic_pointer_cast<GenericObject>(Lobj.value_ptr<T>());
    else
      return boost::shared_ptr<GenericObject>();
  }
  virtual bool set_value(LuabindObject& Lobj,
			 const std::string& Vname, 
			 const boost::shared_ptr<GenericObject>& Obj)
  {
    boost::shared_ptr<T> v = boost::dynamic_pointer_cast<T>(Obj);
    if(v) {
      Lobj.set_value(Vname, v);
      return true;
    }
    return false;
  }
  virtual bool set_index(LuabindObject& Lobj,
			 int Vidx, 
			 const boost::shared_ptr<GenericObject>& Obj)
  {
    boost::shared_ptr<T> v = boost::dynamic_pointer_cast<T>(Obj);
    if(v) {
      Lobj.set_index(Vidx, v);
      return true;
    }
    return false;
  }
};

}
#endif
