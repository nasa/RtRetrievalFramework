from __future__ import absolute_import
from builtins import object
from nose.tools import *
from .xco2_bias_match_matlab_example import *
try:
    import scipy.io as sio
except ImportError:
    pass                        # Ignore trouble importing, we just can't
                                # run unit test
from numpy.testing import *
import os.path
from nose.plugins.skip import Skip, SkipTest
import warnings

class TestXCO2BiasMatchMatlabExample(object):
    def __init__(self):
        if(not have_full_physics_swig):
            return
# Only read data if we are on a system with the data available
        bdir = os.environ.get("L2_HESSIAN_TEST_DIR", 
                 "/groups/algorithm/l2_fp/estimating_bias_and_hessian")
        if(os.path.exists(bdir)):
            self.d = XCO2BiasMatchMatlabExample(bdir +
                                                "/l2_fd_20090605030254.h5")
            self.matlab_result = sio.loadmat(bdir + "/bias_dump.mat")
        else:
            warnings.warn("We don't have the test data available at %s, so we are skipping the XCO2BiasMatchMatlabExample tests." % bdir)
            self.d = None

    def setup(self):
        '''This is run before every test in this class. If we don't have data,
        then we silently skip the test.'''
        if(not have_full_physics_swig):
            raise SkipTest
        if(self.d is None):
            raise SkipTest

    def test_jacobian(self):
        assert_array_almost_equal(self.d.jacobian, self.matlab_result["K"])

    def test_finite_diff_step_size(self):
        assert_array_almost_equal(self.d.finite_diff_step_size, 
                                  self.matlab_result['perturb'][:,0])

    def test_apriori_covariance(self):
        assert_array_almost_equal(self.d.apriori_covariance,
                                  self.matlab_result['Sa'])

    def test_measurement_uncertainty(self):
        assert_array_almost_equal(self.d.measurement_uncertainty,
                                  np.sqrt(np.diag(self.matlab_result['Se'])))

    def test_ma(self):
        assert_array_almost_equal(self.d.ma, 
                                  self.matlab_result['Ma'])

    def test_jacobian_unit_free(self):
        assert_array_almost_equal(self.d.phi0, self.matlab_result['KS'])

    def test_hessian_unit_free(self):
        assert_array_almost_equal(self.d.hessian_unit_free, 
                                  self.matlab_result['HeS'])

    def test_tinv(self):
        assert_array_almost_equal(self.d.tinv, self.matlab_result['Tinv'])

    def test_bias(self):
        # Note that precision of this is a bit less than the default 7 decimal
        # points. I looked at this, and believe this is just a rounding issue.
        # The final xco2_bias does agree to 6 digits, so this seems like it
        # is ok
        assert_array_almost_equal(self.d.bias, self.matlab_result['bias'], 3)

    def test_xco2_bias(self):
        assert_almost_equal(self.d.xco2_bias, self.matlab_result['XCO2_bias'], 
                            6)



    
    
