from nose.tools import *
from full_physics import *
from nose.plugins.skip import Skip, SkipTest
import numpy as np
import numpy.testing as nptest

def test_rayleigh_greek_moment():
    '''Test RayleighGreekMoment'''
    if(not have_full_physics_swig):
        raise SkipTest
    expected = np.array([[1,     0, 0, 0, 0, 0],
                         [1e-11, 0, 0, 1.3968144385817844, 0, 0],
                         [0.47936288771635682, 2.8761773262981407, 0, 0, 
                          1.1741944765321404, 0]])
    nptest.assert_array_almost_equal(RayleighGreekMoment.array(), expected)
    
