from __future__ import absolute_import
from __future__ import division
from builtins import range
from builtins import object
from past.utils import old_div
import numpy as np
from numpy.linalg import inv
from abc import *
from .cached_prop import *

class XCO2Bias(object):
    '''This class calculates an estimate of the XCO2 bias for the Level 2
    full physics calculation of XCO2.

    This calculation is described in the document 
    doc/Estimating_bias_and_Hessian.pdf in the Level 2 source code. This was
    original implemented in matlab. Because of the size of the test data,
    we don't have this checked into the Level 2 source code. But it is
    available in the directory 
    /groups/algorithm/l2_fp/estimating_bias_and_hessian.

    This class is a base class that leaves out details of how exactly 
    the perturbed Jacobians etc. are generated. You will want to use
    one of the classes derived from this one, rather than using XCO2Bias
    directly. In particular, we duplicate the example used in the original
    matlab code in the class XCO2BiasMatchMatlabExample, and we have a
    real calculation for a Level 2 run using XCO2BiasL2Run.
    '''

    @abstractproperty
    def xco2_averaging_kernel(self):
        '''Return the XCO2 averaging kernel'''
        pass

    @abstractproperty
    def xco2_portion_of_bias(self):
        '''In general, the xco2 averaging kernel only applies to the subset
        of the state vector that has the XCO2 values. This should return
        the subset of self.bias that the xco2 averaging kernel applies to.'''
        pass

    @abstractproperty
    def jacobian(self):
        '''Return the unperturbed jacobian'''
        pass

    @abstractproperty
    def jacobian_perturbed(self):
        '''Return the jacobian, for each of the state vector perturbations.
        This is a jac_ncol x jac_nrow x jac_ncol, where [0,:, :] is the
        jacobian for the first state vector element perturbed, [1, :, :]
        is the second, and so on.'''
        pass

    @abstractproperty
    def finite_diff_step_size(self):
        '''Return the step size to use in the finite difference. This is called
        delta_k in the Estimating_bias_and_Hessian.pdf document.'''
        pass

    @abstractproperty
    def apriori_covariance(self):
        '''Return the a priori covariance matrix.'''
        pass

    @abstractproperty
    def measurement_uncertainty(self):
        '''Return measurement uncertainty. This is called sigma_epsilon in
        the Estimating_bias_and_Hessian.pdf document.'''
        pass

    @cached_property
    def se(self):
        return self.measurement_uncertainty * self.measurement_uncertainty

    @cached_property
    def gain(self):
        sa = self.apriori_covariance
        k = self.jacobian
        se = self.se
        return np.dot(inv(inv(sa) + np.dot(old_div(np.transpose(k), se), k)), 
                      old_div(np.transpose(k), se))

    @cached_property
    def avg_kernel(self):
        return np.dot(self.gain, self.jacobian)

    @cached_property
    def mspe(self):
        g = self.gain
        se = self.se
        ak = self.avg_kernel
        sa = self.apriori_covariance
        t = ak - np.eye(ak.shape[0])
        return (np.dot(np.dot(t, sa), np.transpose(t)) +
                np.dot(g * se, np.transpose(g)))

    @cached_property
    def sigma_alpha(self):
        '''This is the a priori uncertainty'''
        return np.sqrt(np.diag(self.apriori_covariance))

    @cached_property
    def ma(self):
        '''This is an intermediate value used in the bias calculation.

        Note that while we usually use np.array, for this set of
        calculations we use np.matrix'''
        r = np.empty_like(self.apriori_covariance)
        for i in range(r.shape[0]):
            r[i,:] = (old_div(self.apriori_covariance[i,:], 
                      (self.sigma_alpha[i] * self.sigma_alpha)))
        return np.mat(r)

    @cached_property
    def tinv(self):
        '''This is an intermediate value used in the bias calculation. This
        is called upsilon on the readme. '''
        p = np.mat(self.phi0)
        return inv(inv(self.ma) + p.T * p)
        
    @cached_property
    def phi0(self):
        '''This is Phi_0 in the paper, which is just a unit free version
        of the unperturbed jacobian.'''
        res = np.empty_like(self.jacobian)
        for i in range(self.jacobian.shape[1]):
            res[:,i] = self.jacobian[:,i] / self.measurement_uncertainty * \
                self.sigma_alpha[i]
        return res

    @cached_property
    def phi_k(self):
        '''This is Phi_k in the paper, which is just a unit free version
        of the perturbed jacobian.'''
        pjac = self.jacobian_perturbed
        res = np.empty_like(pjac)
        for i in range(pjac.shape[2]):
            res[:, :,i] = pjac[:, :,i] / self.measurement_uncertainty * \
                self.sigma_alpha[i]
        return res

    @cached_property
    def hessian_unit_free(self):
        res = np.empty((self.jacobian.shape[0], self.jacobian.shape[1],
                        self.jacobian.shape[1]))
        delta = self.finite_diff_step_size
        for j in range(res.shape[1]):
            for k in range(res.shape[2]):
                res[:,j,k] = (((self.phi_k[k, :, j] - self.phi0[:,j]) /
                               (old_div(delta[k], self.sigma_alpha[k])) *
                               (old_div(delta[k], (delta[j] + delta[k])))
                            ) +
                            ((self.phi_k[j, :, k] - self.phi0[:,k]) /
                             (old_div(delta[j], self.sigma_alpha[j])) *
                             (old_div(delta[j], (delta[j] + delta[k])))
                             ))
        return res

    @cached_property
    def bias(self):
        h = self.hessian_unit_free
        m = self.ma
        t = self.tinv
        res = np.zeros((self.jacobian.shape[1], 1))
        siga = np.mat(np.diag(self.sigma_alpha))
        I = np.mat(np.eye(res.shape[0]))
        p = np.mat(self.phi0)
        for j in range(res.shape[0]):
            hj = np.mat(h[:,:,j])
            hph = hj.T * p + p.T * hj
            temp1 = siga * t * hph * (I - t * p.T * p) * \
                (2 * m * p.T * p + I) * t * siga
            temp2 = siga * t * (-hph * t * p.T + hj.T) * p * t * siga
            tv = temp1[:,j] / 2 / siga[j,j] + temp2[:,j] / 2 / siga[j,j]
            res = res + tv
        return res

    @cached_property
    def xco2_bias(self):
        '''Return the estimate of the XCO2 bias, in parts per million (ppm)'''
        return np.dot(self.xco2_averaging_kernel, 
                      self.xco2_portion_of_bias)[0,0] * 1e6
        
        
