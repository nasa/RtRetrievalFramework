from __future__ import division
# We don't normally need to test all the python interface to our C++ classes.
# We already test the C++ classes in our C++ unit tests, and the interface
# is automatically generated.
#
# However we have these tests in place to check the underlying functionality
# of SWIG (so for example, test that directors are handled correctly.

from nose.tools import *
from full_physics import *
from numpy.testing import *
import numpy as np
from nose.plugins.skip import Skip, SkipTest
import os

test_data = os.path.dirname(__file__) + "/../../unit_test_data/"

# Start some tests using HeritageFile. This is a particularly simple class,
# because it doesn't depend on anything else.

def test_swig_basic_call():
    '''Simple test that we can call C++ through python and get back
    a results.'''
    if(not have_full_physics_swig):
        raise SkipTest
    f = HeritageFile(test_data + "heritage_file_test.run")
    assert f.value_int("ALGORITHMS/points_sun") == 10000

def test_doxygen_doc():
    '''Test that the doxygen documentation goes through'''
    if(not have_full_physics_swig):
        raise SkipTest
    f = HeritageFile(test_data + "heritage_file_test.run")
    assert f.__doc__.find(
'''    This class reads the heritage file formats.

    We read both the configuration file and the matrix files (there are
    similar enough in format that it makes sense to combine these two).''') > 0

@raises(RuntimeError)
def test_exception_handled():
    '''Test that a C++ exception gets translated to a RuntimeError'''
    if(not have_full_physics_swig):
        raise SkipTest
    f = HeritageFile("/home/smyth/Local/Level2/unit_test_data/heritage_file_test.run")
    f.value_bool("WINDOW_INFO/spectral_window_file")

def my_func(self):
    return self.value_int("ALGORITHMS/points_sun") + 10

if(have_full_physics_swig):
    HeritageFile.my_func = my_func

def test_extend_class_in_python():
    '''Test that we can add pure python functions to a class.'''
    if(not have_full_physics_swig):
        raise SkipTest
    f = HeritageFile(test_data + "heritage_file_test.run")
    assert f.my_func() == 10000 + 10

def test_returning_numpy_array():
    '''Test returning a numpy array'''
    if(not have_full_physics_swig):
        raise SkipTest
    f = HeritageFile(test_data + "old_ascii/solar_cont_v1.dat")
    assert_array_equal(f.data, 
                       [[  1.00000000e+00,   8.83596000e+21],
                        [  2.00000000e+00, -9.48206000e+20],
                        [  3.00000000e+00,  -1.51700000e+22],
                        [  4.00000000e+00,   1.74114000e+22],
                        [  5.00000000e+00, -7.73485000e+21],
                        [  6.00000000e+00,   1.23130000e+21]])

def test_second_class():
    '''Test is a second class. This makes sure the module initialization etc.
    get handled correctly for multiple classes.'''
    if(not have_full_physics_swig):
        raise SkipTest
    assert conversion(Unit("m"), Unit("cm")) == 100.0
    u = Unit("m")
    assert_almost_equal(u.conversion_to_si, 1.0)
    assert u.name == "m"
    # Not working yet
    #print u.base_unit_powers

def test_header_only_class():
    '''Test handling of a class with a .h but not a .cc file. This is 
    mainly a test of the Makefile, SWIG doesn't really act any differently
    for .h vs .h and .cc files'''
    if(not have_full_physics_swig):
        raise SkipTest
    t = Exception("My exception")
    assert t.what() == "My exception"

def test_class_using_another():
    '''Test that import is working properly when a class uses another one. 
    We also test that the python operators (e.g., __mul__) get handled
    correctly.'''
    if(not have_full_physics_swig):
        raise SkipTest
    d = DoubleWithUnit(10, "m")
    d2 = d.convert("cm")
    assert_almost_equal(d2.value, 10 * 100.0)
    d2.value = 4
    d2.units = Unit("m")
    assert_almost_equal(d2.convert("m").value, 4.0)
    d2 = d2.convert("cm")
    d3 = d + d2
    assert_almost_equal(d3.convert("m").value, 14.0)
    d3 = d - d2
    assert_almost_equal(d3.convert("m").value, 6.0)
    d3 = d / d2
    assert_almost_equal(d3.convert("dimensionless").value, 10.0 / 4)
    d3 = d * d2
    assert_almost_equal(d3.convert("m^2").value, 40.0)

def test_passing_numpy_array():
    '''Test that passes a numpy array to a function'''
    if(not have_full_physics_swig):
        raise SkipTest
    a = np.array([[1, 2],[3,4],[5,6]])
    t = ArrayWithUnit_double_2(a, "m")
    assert_almost_equal(t.value, a)
    b = np.array([[10, 20],[3,4],[55,66]])
    t.value = b
    assert_almost_equal(t.value, b)
    t = ArrayWithUnit_double_2(a[0:2,:], "m")
    assert_almost_equal(t.value, [[1,2],[3,4]])
    t = ArrayWithUnit_double_2(a[0::2,:], "m")
    assert_almost_equal(t.value, [[1,2],[5,6]])
    t = ArrayWithUnit_double_1(a[0:2,1], "m")
    assert_almost_equal(t.value, [2,4])
