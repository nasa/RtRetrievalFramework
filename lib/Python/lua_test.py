from builtins import object
from nose.tools import *
from full_physics import *
from nose.plugins.skip import Skip, SkipTest
import os

test_data = os.path.dirname(__file__) + "/../../unit_test_data/"

if(not have_full_physics_swig):
    class LuaCallback(object):
        pass

# A callback that takes no arguments
class MyCallback0(LuaCallback):
    def __init__(self, ls):
        LuaCallback.__init__(self, ls)
        self.ls = ls
        
    def call(self, obj1,obj2,obj3,obj4,obj5,obj6,obj7,obj8,obj9,obj10):
        return LuabindObject(self.ls, 1)
        
# A callback that takes one arguments
class MyCallback1(LuaCallback):
    def __init__(self, ls):
        LuaCallback.__init__(self, ls)
        self.ls = ls
                
    def call(self, obj1,obj2,obj3,obj4,obj5,obj6,obj7,obj8,obj9,obj10):
        obj1.v3 = "hi there"
        return LuabindObject(self.ls, 2)

class TestLua(object):
    def setUp(self):
        if(not have_full_physics_swig):
            raise SkipTest
        self.ls = LuaState(test_data)

    # Test basic reading and writing of object to Lua
    def test_read_and_write(self):
        self.ls.run("test_var = {}\n")
        test_var = self.ls.globals.test_var
        test_var.v1 = "hi there"
        test_var.v2 = 1
        test_var.v3 = 2.0
        t = HeritageFile(test_data + "heritage_file_test.run")
        assert t.value_int("ALGORITHMS/points_sun") == 10000
        test_var.test_obj = t
        assert self.ls.globals.test_var.v1 == "hi there"
        assert self.ls.globals.test_var.v2 == 1
        assert self.ls.globals.test_var.v3 == 2.0
        assert self.ls.globals.test_var.test_obj.value_int("ALGORITHMS/points_sun") == 10000

        test_list = ["a", "b", "abc"]
        test_var.test_list = test_list
        for l_idx, l_val in enumerate(test_list):
            assert self.ls.globals.test_var.test_list[l_idx+1] == test_list[l_idx]

        test_dict = {"k1": "abc", "k2": "cde", 1: "blah"}
        test_var.test_dict = test_dict
        for l_key, l_val in list(test_dict.items()):
            assert self.ls.globals.test_var.test_dict[l_key] == test_dict[l_key]

    # Test reading a config file
    def test_config(self):
        self.ls.do_file(test_data + "config.lua")
        assert self.ls.globals.config.atmosphere.pressure.surface_pressure.value.value == 96716.6249

    # Test a callback function object. Note you don't normally use this 
    # directly, but rather use the callback tested in the next function
    def test_callback_object(self):
        self.ls.run("test_var = {}\n")
        test_var = self.ls.globals.test_var
        test_var.f0 = MyCallback0(self.ls)
        test_var.f1 = MyCallback1(self.ls)
        global test_v
        test_v = 0
        assert test_v == 0
        self.ls.run("test_v = test_var.f0()")
        assert self.ls.globals.test_v == 1
        self.ls.run("test_v = test_var:f1()")
        assert self.ls.globals.test_v == 2
        assert test_var.v3 == "hi there"

    def f0(self):
        global test_v
        test_v = 1

    def f1(self, v):
        global test_v
        test_v = v

    def f2(self, v1, v2):
        global test_v
        test_v = v1 + v2

    def f3(self, v1, v2, v3):
        global test_v
        test_v = v1 + v2 + v3

    def f3_return(self, v1, v2, v3):
        return v1 + v2 + v3

    # Test passing any callable object to Lua, to make sure the callbacks work
    def test_callback(self):
        g = self.ls.globals
        g.f0 = self.f0
        g.f1 = self.f1
        g.f2 = self.f2
        g.f3 = self.f3
        g.f3_return = self.f3_return
        global test_v
        test_v = 0
        assert test_v == 0
        self.ls.run("f0(nil,nil,nil,nil,nil,nil,nil,nil,nil)")
        assert test_v == 1
        self.ls.run("f1(3.5)")
        assert test_v == 3.5
        self.ls.run("f2(3, 4)")
        assert test_v == 3 + 4
        self.ls.run("f3(3, 4, 5)")
        assert test_v == 3 + 4 + 5
        self.ls.run("val = f3_return(3, 4, 5)")
        assert g.val == 3 + 4 + 5

    # Test calling a Lua function.
    def test_luafunc(self):
        self.ls.run('''
                    function test_func(v1, v2)
                       return v1 + v2 * 2
                    end''')
        tf = self.ls.globals.test_func
        assert tf(1, 2) == 1 + 2 * 2
        self.ls.run('''
                    function test_func2()
                       return "blah"
                    end''')
        tf2 = self.ls.globals.test_func2
        assert tf2() == "blah"
        
