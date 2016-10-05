from nose.tools import *
from full_physics import *
from nose.plugins.skip import Skip, SkipTest
import os

test_data = os.path.dirname(__file__) + "/../../unit_test_data/"

def test_heritage_file():
    if(not have_full_physics_swig):
        raise SkipTest
    f = HeritageFile(test_data + "heritage_file_test.run")
    assert f.value_int("ALGORITHMS/points_sun") == 10000
    
