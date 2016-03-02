from future import standard_library
standard_library.install_aliases()
# Test that we set up pickling correctly in the SWIG code.

from nose.tools import *
from full_physics import *
from nose.plugins.skip import Skip, SkipTest
import os
import pickle

test_data = os.path.dirname(__file__) + "/../../unit_test_data/"

def test_hdf_file_pickle():
    '''Test pickling of HdfFile class.'''
    if(not have_full_physics_swig):
        raise SkipTest
    f = HdfFile(test_data + "l1b.h5")
    t = pickle.dumps(f, pickle.HIGHEST_PROTOCOL)
    f2 = pickle.loads(t)
    assert f.file_name == f2.file_name
