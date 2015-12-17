#include "unit_test_support.h"
#include <iostream>
extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}
#include <luabind/luabind.hpp>
#include "level_1b.h"
#include "register_lua.h"
#include "hdf_file.h"
using namespace FullPhysics;

namespace FullPhysics {
  void fake_func();
}
class LuaFixture: public GlobalFixture {
public:
  LuaFixture() 
  {ls = luaL_newstate(); luaL_openlibs(ls); luabind::open(ls);}
  virtual ~LuaFixture() {   lua_close(ls); }
  lua_State *ls;
};

BOOST_FIXTURE_TEST_SUITE(lua, LuaFixture)

BOOST_AUTO_TEST_CASE(hello_world)
{
  return;			// Don't normally run this, since it
				// doesn't produce any kind of test
				// results. But leave in place for
				// basic diagnose of Lua problems.
  // Basic hello world, so we make sure we have the basic pieces in
  // place to call Lua.
  std::string fname = test_data_dir() + "lua/hello_world.lua";
  if(luaL_loadfile(ls, fname.c_str()) ||
     lua_pcall(ls, 0, 0, 0)) {
    // If we are here, then there is an error message on the stack
    std::cerr << "Lua error: " << lua_tostring(ls, -1) << "\n";
    lua_pop(ls, 1);		// Remove error message from stack.
  }
}

BOOST_AUTO_TEST_CASE(luabind_simple_func)
{
  // Create a simple add function on Lua, and use Luabind to call.

  luaL_dostring(ls,
		"function add(first, second)\n"
		"  return first + second\n"
		"end\n"
		);
  BOOST_CHECK_EQUAL((int) luabind::call_function<int>(ls, "add", 2, 3),
  		    5);
}

int my_broken_add(int x, int y) { return x + y + 2; }
BOOST_AUTO_TEST_CASE(luabind_simple_func_to_lua)
{
  // Now we define a function in C++, pass it to Lua, and have
  // Lua call it.
  luabind::module(ls)[luabind::def("my_broken_add", my_broken_add)];
  luaL_dostring(ls,
		"function add(first, second)\n"
		"  return my_broken_add(first, second)\n"
		"end\n"
		);
  BOOST_CHECK_EQUAL((int) luabind::call_function<int>(ls, "add", 2, 3),
  		    5 + 2);
}

class TestClass {
public:
  TestClass(int x) : x_(x) { }
  virtual ~TestClass() {}
  double y;
  virtual double add(double z) { return x_ + y + z; }
private:
  int x_;
};

BOOST_AUTO_TEST_CASE(luabind_simple_class)
{
  // Wrap a class, and use in in Lua.
  luabind::module(ls)
    [luabind::class_<TestClass>("TestClass")
     .def(luabind::constructor<int>())
     .def("add", &TestClass::add)
     .def_readwrite("y", &TestClass::y)];
  luaL_dostring(ls,
		"function add(first, second, third)\n"
                "  t = TestClass(first)\n"
                "  t.y = second\n"
                "  return t:add(third)\n"
		"end\n"
		);
  BOOST_CHECK_CLOSE((double) luabind::call_function<double>(ls, "add", 
							    1, 2.5, 3.5),
  		    1 + 2.5 + 3.5, 1e-8);
}

BOOST_AUTO_TEST_CASE(luabind_return_class)
{
  // Wrap a class, and use in in Lua. Have Lua return an instance of the
  // class that it creates.
  luabind::module(ls)
    [luabind::class_<TestClass,boost::shared_ptr<TestClass> >("TestClass")
     .def(luabind::constructor<int>())
     .def("add", &TestClass::add)
     .def_readwrite("y", &TestClass::y)];
  luaL_dostring(ls,
		"function create_class(first)\n"
                "  t = TestClass(first)\n"
                "  return t\n"
		"end\n"
		);
  boost::shared_ptr<TestClass> t(luabind::call_function<boost::shared_ptr<TestClass> >(ls, "create_class", 1));
  t->y = 2.5;
  BOOST_CHECK_CLOSE(t->add(3.5), 1 + 2.5 + 3.5, 1e-8);
}

class TestClassDerived : public TestClass {
public:
  TestClassDerived(int x, double x2) : TestClass(x), x2_(x2) { }
  virtual ~TestClassDerived() {}
  virtual double add(double z) { return TestClass::add(z) + x2_; }
private:
  double x2_;
};

BOOST_AUTO_TEST_CASE(luabind_return_derived_class)
{
  // Wrap a class, and use in in Lua. Have Lua return an instance of the
  // class that it creates.
  luabind::module(ls)
    [luabind::class_<TestClass,boost::shared_ptr<TestClass> >("TestClass")
     .def(luabind::constructor<int>())
     .def("add", &TestClass::add)
     .def_readwrite("y", &TestClass::y)];
  luabind::module(ls)
    [
     luabind::class_<TestClassDerived, TestClass, 
              boost::shared_ptr<TestClass> >("TestClassDerived")
     .def(luabind::constructor<int, double>())
     .def("add", &TestClassDerived::add)
     ];
  luaL_dostring(ls,
		"function create_class(first)\n"
                "  t = TestClassDerived(first, 11.5)\n"
                "  return t\n"
		"end\n"
		);
  boost::shared_ptr<TestClass> t(luabind::call_function<boost::shared_ptr<TestClass> >(ls, "create_class", 1));
  t->y = 2.5;
  BOOST_CHECK_CLOSE(t->add(3.5), 1 + 2.5 + 3.5 + 11.5, 1e-8);
}

BOOST_AUTO_TEST_CASE(register_class)
{
  RegisterLua::register_lua(ls);
  std::string s = test_data_dir() + "l1b.h5";
  std::string sid = "20090725020316";
  std::ostringstream os;
  os << "s = '" << s << "'\n"
     << "l1b_hdf = HdfFile(s)"
     << "sid = AcosSoundingId(l1b_hdf, '" << sid << "', AcosSoundingId.P_SOUNDING)\n"
     << "l1b = Level1bAcos(l1b_hdf, sid)\n";
  luaL_dostring(ls, os.str().c_str());
  boost::shared_ptr<Level1b> l1b =
    luabind::object_cast<boost::shared_ptr<Level1b> >
    (luabind::globals(ls)["l1b"]);
  BOOST_CHECK_EQUAL(l1b->number_spectrometer(), 3);
  BOOST_CHECK_CLOSE(l1b->radiance(0).data()(402 + 10), 
		    5.0381521532472107e-07, 1e-4);
}

class FuncWrap {
public:
  FuncWrap() {}
  void operator()(const luabind::object& F) {  }
};

BOOST_AUTO_TEST_CASE(set_value)
{
  RegisterLua::register_lua(ls);
  luabind::module(ls)[luabind::class_<FuncWrap, boost::shared_ptr<FuncWrap> >("FuncWrap").def(luabind::constructor<>()).def("__call", &FuncWrap::operator())];
  std::ostringstream os;
  os << "s = { }\n";
  luaL_dostring(ls, os.str().c_str());
  std::string s = test_data_dir() + "l1b.h5";
  boost::shared_ptr<HdfFile> h(new HdfFile(s));
  luabind::globals(ls)["s"]["t"] = h;
  luaL_dostring(ls, "s.t2 = s.t\n");
  boost::shared_ptr<HdfFile> h2 =
    luabind::object_cast<boost::shared_ptr<HdfFile> >(luabind::globals(ls)["s"]["t2"]);
  BOOST_CHECK_EQUAL(h2->file_name(), h->file_name());
  boost::shared_ptr<FuncWrap> f(new FuncWrap);
  luabind::globals(ls)["s"]["f"] = f;
  luaL_dostring(ls, "s:f()\n");
}

BOOST_AUTO_TEST_SUITE_END()
