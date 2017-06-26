from __future__ import absolute_import
from .xco2_bias import *
from full_physics import *

class XCO2BiasMatchMatlabExample(XCO2Bias):
    '''This is an instance of XCO2 that reads all the Jacobians etc. from
    an existing Level 2 like HDF file. This is used to duplicate the example
    supplied with the original matlab code. Once we have established that
    everything is working, this class is not very useful - this is really
    only meant to support testing.
    '''
    def __init__(self, fname = "/groups/algorithm/l2_fp/estimating_bias_and_hessian/l2_fd_20090605030254.h5"):
        self.hfile = HdfFile(fname)
        self.jac = self.hfile.read_double_3d("/RetrievalResults/jacobian")

    @property
    def jacobian(self):
        '''Return the unperturbed jacobian'''
        return self.jac[0,:,:]

    @property
    def jacobian_perturbed(self):
        '''Return the jacobian, for each of the state vector perturbations.
        This is a jac_ncol x jac_nrow x jac_ncol, where [0,:, :] is the
        jacobian for the first state vector element perturbed, [1, :, :]
        is the second, and so on.'''
        return self.jac[1:,:,:]

    @cached_property
    def finite_diff_step_size(self):
        '''Return the step size to use in the finite difference. This is called
        delta_k in the Estimating_bias_and_Hessian.pdf document.'''
        # Just hardcode this for our test. This is what was done in the
        # matlab code.
        delta = [1e-7] * 20
        delta.extend([0.01, 10.0])
        delta.extend([0.1] * (20 * 4 + 1))
        delta.extend([0.001, 1e-6, 1e-3, 1e-6,1e-3, 1e-6])
        delta.extend([5e-4] * 3)
        return delta

    @cached_property
    def apriori_covariance(self):
        '''Return the a priori covariance matrix.'''
        return self.hfile.read_double_3d("/RetrievalResults/apriori_covariance_matrix")[0,:,:]

    @cached_property
    def measurement_uncertainty(self):
        '''Return measurement uncertainty. This is called sigma_epsilon in
        the Estimating_bias_and_Hessian.pdf document.'''
        return self.hfile.read_double_2d("/SpectralParameters/measured_radiance_uncert")[0,:]

    @cached_property
    def xco2_averaging_kernel(self):
        return self.hfile.read_double_2d("/RetrievalResults/xco2_pressure_weighting_function")[0,:]

    @property
    def xco2_portion_of_bias(self):
        '''In general, the xco2 averaging kernel only applies to the subset
        of the state vector that has the XCO2 values. This should return
        the subset of self.bias that the xco2 averaging kernel applies to.'''
        return self.bias[0:self.xco2_averaging_kernel.shape[0], 0]




