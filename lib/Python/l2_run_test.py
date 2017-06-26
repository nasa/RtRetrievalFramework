from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from nose.tools import *
from .l2_run import *
from nose.plugins.skip import Skip, SkipTest
import os
import pickle

gosat_config = os.path.dirname(__file__) + "/../../input/gosat/config/config.lua"
tccon_small_set = os.path.dirname(__file__) + "/../../test/tccon_small_set/"
met_file = tccon_small_set + "acos_EcmB2900_tccon_5_good_qual.h5"
spectrum_file = tccon_small_set + "acos_L1bB2900_tccon_5_good_qual.h5"
sounding_id = "20100223034944"

def test_l2_run():
    '''Basic test of L2Run object.'''
    if(not have_full_physics_swig):
        raise SkipTest
    r = L2Run(gosat_config, sounding_id, met_file, spectrum_file)

def test_l2_run_pickle():
    if(not have_full_physics_swig):
        raise SkipTest
    r = L2Run(gosat_config, sounding_id, met_file, spectrum_file)
    t = pickle.dumps(r, pickle.HIGHEST_PROTOCOL)
    r2 = pickle.loads(t)


