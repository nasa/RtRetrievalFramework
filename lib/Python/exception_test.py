from nose.tools import *
from full_physics import *
from nose.plugins.skip import Skip, SkipTest

def test_exception():
    if(not have_full_physics_swig):
        raise SkipTest
    t = Exception("test")
    assert t.what() == "test"
        
