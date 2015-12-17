#include "unit_test_support.h"
#include "lua_state.h"
#include "lua_callback.h"
#include "hdf_file.h"
using namespace FullPhysics;

class LuaCallbackTest: public LuaCallback {
public:
  LuaCallbackTest(int Val, const LuaState& Ls) 
    : LuaCallback(Ls), val(Val), ls(new LuaState(Ls)) {}
  virtual ~LuaCallbackTest() {}
  virtual boost::shared_ptr<LuabindObject> 
  call(const boost::shared_ptr<LuabindObject>& Obj,
       const boost::shared_ptr<LuabindObject>& Obj2,
       const boost::shared_ptr<LuabindObject>& Obj3,
       const boost::shared_ptr<LuabindObject>& Obj4,
       const boost::shared_ptr<LuabindObject>& Obj5,
       const boost::shared_ptr<LuabindObject>& Obj6,
       const boost::shared_ptr<LuabindObject>& Obj7,
       const boost::shared_ptr<LuabindObject>& Obj8,
       const boost::shared_ptr<LuabindObject>& Obj9,
       const boost::shared_ptr<LuabindObject>& Obj10)
  { Obj->set_value("call_back", val); return boost::shared_ptr<LuabindObject>(new LuabindObject(LuabindObject::nil(ls))); }
private:
  int val;
  boost::shared_ptr<LuaState> ls;
};
BOOST_FIXTURE_TEST_SUITE(lua_state, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<LuaState> ls =
    LuaState::load_file(test_data_dir() + "lua_file_test.lua");
  LuabindObject lf = ls->globals();
  BOOST_CHECK_EQUAL(lf["l1b_hdf"].value_ptr<HdfFile>()->file_name(), 
		    lf["l1b_fname"].value<std::string>());
  boost::shared_ptr<GenericObject> gobj = 
    lf["l1b_hdf"].value_generic_object();
  BOOST_CHECK(gobj);
  lf.set_value("t", "hi there");
  lf.set_value("t2", 10);
  lf.set_value("t3", 10.5);
  gobj = lf["t"].value_generic_object();
  BOOST_CHECK(!gobj);
  BOOST_CHECK_EQUAL(lf["t"].value<std::string>(), "hi there");
  BOOST_CHECK_EQUAL(lf["t2"].value<int>(), 10);
  BOOST_CHECK_CLOSE(lf["t3"].value<double>(), 10.5, 1e-8);
  BOOST_CHECK(lf["t"].is_string());
  BOOST_CHECK(lf["t2"].is_number());
  BOOST_CHECK(lf["t3"].is_number());
  BOOST_CHECK(lf["t4"].is_nil());
  ls->run("v = {}");
  LuabindObject v = lf["v"];
  boost::shared_ptr<LuaCallback> lc(new LuaCallbackTest(34, *ls));
  v.set_value("func", lc);
  ls->run("v:func()");
  BOOST_CHECK_EQUAL(v["call_back"].value<int>(), 34);
  ls->run("function test_func(v1, v2)\n  return v1 + 2 * v2\nend\n");
  BOOST_CHECK_EQUAL(lf["test_func"].call(1, 2).value<int>(),
		    1 + 2 * 2);
}

BOOST_AUTO_TEST_CASE(array)
{
  boost::shared_ptr<LuaState> ls =
    LuaState::load_file(test_data_dir() + "lua_file_test.lua");
  LuabindObject lf = ls->globals();
  blitz::Array<double, 1> t(3);
  t = 1, 2, 3;
  lf.set_value("t", t);
  blitz::Array<double, 1> tres(lf["t"].value<blitz::Array<double, 1> >());
  BOOST_CHECK_MATRIX_CLOSE(tres, t);
}

BOOST_AUTO_TEST_SUITE_END()
