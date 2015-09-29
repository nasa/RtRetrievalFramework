#include "luabind_object.h"
#include <boost/foreach.hpp>
using namespace FullPhysics;

std::vector<boost::shared_ptr<BridgeLuabindAndGenericBase> > 
LuabindObject::base_bridge;

std::vector<boost::shared_ptr<BridgeLuabindAndGenericBase> > 
LuabindObject::derived_bridge;

//-----------------------------------------------------------------------
/// Return value as a GenericObject, or return a null pointer if we
/// can't convert to a GenericObject. This is useful for use with
/// Python. 
//-----------------------------------------------------------------------

boost::shared_ptr<GenericObject> LuabindObject::value_generic_object() const
{
  boost::shared_ptr<GenericObject> res;
  BOOST_FOREACH(const boost::shared_ptr<BridgeLuabindAndGenericBase>& v,
		base_bridge) {
    res = v->value(*this);
    if(res)
      return res;
  }
  return res;
}

//-----------------------------------------------------------------------
/// Set value as a GenericObject, returning true if this succeeds. This 
/// is useful for use with Python. 
//-----------------------------------------------------------------------

void LuabindObject::set_value
(const std::string& Vname, 
 const boost::shared_ptr<GenericObject>& V)
{
  BOOST_FOREACH(const boost::shared_ptr<BridgeLuabindAndGenericBase>& v,
		derived_bridge) {
    if(v->set_value(*this, Vname, V))
      return;
  }
  BOOST_FOREACH(const boost::shared_ptr<BridgeLuabindAndGenericBase>& v,
		base_bridge) {
    if(v->set_value(*this, Vname, V))
      return;
  }
  throw Exception("Unrecognized type in LuabindObject::set_value"); 
}

//-----------------------------------------------------------------------
/// Set value as a GenericObject, returning true if this succeeds. This 
/// is useful for use with Python. 
//-----------------------------------------------------------------------

void LuabindObject::set_index
(int Vidx, 
 const boost::shared_ptr<GenericObject>& V)
{
  BOOST_FOREACH(const boost::shared_ptr<BridgeLuabindAndGenericBase>& v,
		derived_bridge) {
    if(v->set_index(*this, Vidx, V))
      return;
  }
  BOOST_FOREACH(const boost::shared_ptr<BridgeLuabindAndGenericBase>& v,
		base_bridge) {
    if(v->set_index(*this, Vidx, V))
      return;
  }
  throw Exception("Unrecognized type in LuabindObject::set_index"); 
}


//-----------------------------------------------------------------------
/// Create a LuabindObject from a GenericObject. This is really meant
/// for use in python, in C++ you can just directly use the templated
/// constructor with the actually type.
//-----------------------------------------------------------------------

boost::shared_ptr<LuabindObject> 
LuabindObject::create_luabind_object
(const boost::shared_ptr<LuaState>& Ls,
 const boost::shared_ptr<GenericObject>& Obj)
{
  boost::shared_ptr<LuabindObject> res;
  BOOST_FOREACH(const boost::shared_ptr<BridgeLuabindAndGenericBase>& v,
		derived_bridge) {
    res = v->luabind_object(Ls, Obj);
    if(res)
      return res;
  }
  BOOST_FOREACH(const boost::shared_ptr<BridgeLuabindAndGenericBase>& v,
		base_bridge) {
    res = v->luabind_object(Ls, Obj);
    if(res)
      return res;
  }
  throw Exception("Unrecognized type in LuabindObject::set_index"); 
}
