#include "array_with_unit.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"

DoubleWithUnit awu_double_1d_read(const ArrayWithUnit<double, 1>& V, int i)
{ return V(i); }
DoubleWithUnit awu_double_2d_read(const ArrayWithUnit<double, 1>& V, int i, int j)
{ return V(i, j); }
DoubleWithUnit awu_double_3d_read(const ArrayWithUnit<double, 1>& V, int i, int j, int k)
{ return V(i, j, k); }

template <class A>
std::string array_with_unit_tostring(const A& V)
{
  std::ostringstream os;
  os << V;
  return os.str();
}

template <class A>
std::string array_with_unit_unit_get(const A& V)
{
  return V.units.name();
}

template <class A>
void array_with_unit_unit_set(A& V, std::string& Unit_name)
{
  V.units = Unit(Unit_name);
}


typedef ArrayWithUnit<double, 1> awud1;
typedef ArrayWithUnit<double, 2> awud2;
typedef ArrayWithUnit<double, 3> awud3;
REGISTER_LUA_CLASS_NAME(awud1, ArrayWithUnit_1d)
.def(luabind::constructor<const blitz::Array<double, 1>&, 
			  const std::string&>())
.def_readwrite("value", &awud1::value)
.property("units", 
          &array_with_unit_unit_get<awud1>,
          &array_with_unit_unit_set<awud1>)
.def("__call", &awu_double_1d_read)
.def("__tostring", &array_with_unit_tostring<awud1>)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(awud2, ArrayWithUnit_2d)
.def(luabind::constructor<const blitz::Array<double, 2>&, 
			  const std::string&>())
.def_readwrite("value", &awud2::value)
.property("units", 
          &array_with_unit_unit_get<awud2>,
          &array_with_unit_unit_set<awud2>)
.def("__call", &awu_double_2d_read)
.def("__tostring", &array_with_unit_tostring<awud2>)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(awud3, ArrayWithUnit_3d)
.def(luabind::constructor<const blitz::Array<double, 3>&, 
			  const std::string&>())
.def_readwrite("value", &awud3::value)
.property("units", 
          &array_with_unit_unit_get<awud3>,
          &array_with_unit_unit_set<awud3>)
.def("__call", &awu_double_3d_read)
.def("__tostring", &array_with_unit_tostring<awud3>)
REGISTER_LUA_END()

#endif
