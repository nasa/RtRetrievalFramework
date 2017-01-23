from __future__ import absolute_import
from builtins import str
from nose.tools import *
import os
from .try_swig_load import *
from .acos_file import *
from nose.plugins.skip import Skip, SkipTest
import numpy as np

test_data = os.path.dirname(__file__) + "/../../unit_test_data/"
l2_full_test_data = os.path.dirname(__file__) + "/../../test/l2_full"

# We can't have both h5py and HdfFile work in the same program (never
# have been able to figure out the underlying reason for this). As a work
# around, only run unit tests here if we don't have swig so the HdfFile test
# isn't being run.
def test_gosat_file():
    if(have_full_physics_swig):
        raise SkipTest
    gosat_file = test_data + "l1b.h5"
    gosat_obj = L1B(gosat_file)
    assert gosat_obj.get_sounding_indexes("20090725015225") == (5,slice(2))
    assert gosat_obj.get_sounding_indexes("20090725015331") == (20,slice(2))
    assert gosat_obj.get_sounding_indexes(20090725015924) == (99,slice(2))
    assert gosat_obj.get_sounding_indexes(20090725021407) == (297,slice(2))

    # Set these so we are guaranteed to read from file
    gosat_obj._data_shape_name_dict = {}
    gosat_obj._default_shape_names = None
    
    assert gosat_obj.get_data_shape('/FootprintGeometry/footprint_stokes_coefficients') == ['Exposure', 'Band', 'Polarization', 'StokesCoefficient']

    sounding_id = 20090725015225
    read_latitude = gosat_obj.get_sounding_info('sounding_latitude', sounding_id)
    expt_latitude = np.array( [[ 69.85140991,  69.85140991],
                                  [ 69.85140991,  69.85140991],
                                  [ 69.85140991,  69.85140991]] )
    assert np.all(abs(read_latitude - expt_latitude) < 1e-6)

    read_latitude = gosat_obj.get_sounding_info('sounding_latitude', sounding_id, average='Polarization')
    expt_latitude = numpy.array( [ 69.85140991,  69.85140991, 69.85140991] )
    assert np.all(abs(read_latitude - expt_latitude) < 1e-6)

    read_latitude = gosat_obj.get_sounding_info('sounding_latitude', str(sounding_id) + "P")
    assert np.all(abs(read_latitude - expt_latitude) < 1e-6)

    assert gosat_obj.get_sounding_time(sounding_id)[0] == (2009, 7, 25, 1, 52, 26, 5, 206, -1)
        
    assert  gosat_obj.get_surface_grouping(sounding_id) == 'land'
    assert gosat_obj.get_build_id() == (2, 6, 1)

def test_oco_file():
    if(have_full_physics_swig):
        raise SkipTest
    oco_file = test_data + "oco2_sim_l1b.h5"
    oco_obj = L1B(oco_file)
    assert oco_obj.get_sounding_indexes("2010090912004075") == (0,2)
    assert oco_obj.get_sounding_indexes("2010090912004076") == (0,3)

    sounding_id = 2010090912004075
    read_latitude = oco_obj.get_sounding_info('sounding_latitude', sounding_id)
    expt_latitude = numpy.array( [ -20.02293015, -20.02293015, -20.02293015] )
    assert np.all(abs(read_latitude - expt_latitude) < 1e-6)

    # leap seconds should be accounted here
    assert oco_obj.get_sounding_time(sounding_id)[0] == (2010, 9, 9, 12, 0, 40, 3, 252, -1)


def test_l2_file():
    if(have_full_physics_swig):
        raise SkipTest
    l2_file = l2_full_test_data + "/out.expected.h5"
    l2_obj = L2(l2_file)

    assert l2_obj.get_sounding_indexes("2006091412272205") == (0,)

    assert l2_obj.get_data_shape('/RetrievalResults/xco2_avg_kernel') == ['Retrieval', 'Level']
        
