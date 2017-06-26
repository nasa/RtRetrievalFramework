from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import map
from builtins import range
from .xco2_bias import *
from .l2_run import *
from .orbit_sim import *
import pickle
import re
from functools import partial

def jacobian_calc(self, num, i):
    print("===============================================")
    print("Calculation %d of %d:" % (i + 1, num))
    print("===============================================")
    print("")
    return self.jacobian_perturbed_i(i)

class XCO2BiasFromL2Run(XCO2Bias):
    '''This is an instance of XCO2Bias that calculates everything from an
    existing Level 2 run.

    You can optionally pass in the name of a file containing a pickled 
    distribution covariance to use instead of the Apriori_covariance used in
    the retrieval'''
    def __init__(self, config_name, sounding_id = None,
                 sounding_id_index = 0, lua_config = None, pool = None,
                 dist_cov_file = None):
        self.r = L2Run.create_apriori_state(config_name, sounding_id, 
                                            sounding_id_index, 
                                            lua_config)
        self.x_a = np.copy(self.r.state_vector.state)
        self.pool = pool
        self.dist_cov_file = dist_cov_file

    @property
    def xco2_true(self):
        '''Determine xco2 true. Note that we need the log file, which isn't
        passed in. We assume that the log file is in same directory as
        the simulator met file (it usually is). We can extend this later
        if that because an issue.'''
        logfile = re.sub(r'.hdf', '.log', 
                       re.sub(r'meteorology', 'scene', self.r.met_file))
        lf = OrbitSimLogFile(logfile, self.r.spectrum_file)
        return lf.data[lf.get_sounding_indexes(self.r.sounding_id).Frame, 0,
                       lf.get_column_index("XCO2")]

    @property
    def sounding_id(self):
        return self.r.sounding_id

    @property
    def state_vector_name(self):
        return self.r.state_vector.state_vector_name

    @property
    def x_sol(self):
        return self.r.output_file.read_double_2d("/RetrievedStateVector/state_vector_result")

    @cached_property
    def jacobian(self):
        '''Return the unperturbed jacobian'''
        cf = self.r.cost_function
        residual, se, jac = cf.cost_function(self.x_a)
        return jac

    @cached_property
    def jacobian_perturbed(self):
        '''Return the jacobian, for each of the state vector perturbations.
        This is a jac_ncol x jac_nrow x jac_ncol, where [0,:, :] is the
        jacobian for the first state vector element perturbed, [1, :, :]
        is the second, and so on.'''
        # We may try to do this in parallel at some point, but for now
        # just do this sequentially
        res = np.empty((self.jacobian.shape[1], self.jacobian.shape[0],
                        self.jacobian.shape[1]))
        # Can't pickle a pool object, so go ahead an grab and remove before
        # we possibly pass though a pool. We'll grab this back after we are
        # done
        pool_save = self.pool
        self.pool = None
        try:
            proc = partial(jacobian_calc, self, res.shape[0])
            if(pool_save is None):
                t = list(map(proc, list(range(res.shape[0]))))
            else:
                t = pool_save.map(proc, list(range(res.shape[0])))
            for i in range(res.shape[0]):
                res[i, :, :] = t[i]
        finally:
            self.pool = pool_save
        return res

    def jacobian_perturbed_i(self, i):
        '''Return the jacobian, with the i_th state vector element 
        perturbed by finite_diff_step_size (zero based).'''
        try:
            x = np.copy(self.x_a)
            x[i] += self.finite_diff_step_size[i]
            cf = self.r.cost_function
            residual, se, jac = cf.cost_function(x)
            return jac
        except RuntimeError:
            print("Jacobian_perturbed_id Failed for i = %d" % i)
            raise


    @cached_property
    def finite_diff_step_size(self):
        '''Return the step size to use in the finite difference. This is called
        delta_k in the Estimating_bias_and_Hessian.pdf document.'''
        return np.array([self.finite_diff_step_size_i(i) for i in range(len(self.x_a))])

    def finite_diff_step_size_i(self, i):
        '''Right now, this is hard coded based on the state vector name.
        This is a bit fragile, we may need to come up with something better.'''
        svname =  self.r.state_vector.state_vector_name[i]
        if(re.search(r"CO2 VMR for Press", svname)):
            return 1e-7
        elif(re.search(r"H2O Scaling factor", svname)):
            return 0.01
        elif(re.search(r"Surface Pressure", svname)):
            return 10.0
        elif(re.search(r"Temperature", svname)):
            return 0.1
        elif(re.search(r'Ground Lambertian \S+ Albedo Parm 1', svname)):
            return 0.001
        elif(re.search(r'Ground Lambertian \S+ Albedo Parm 2', svname)):
            return 1e-6
        elif(re.search(r'Instrument Dispersion \S+ Offset', svname)):
            return 5e-4
        elif(re.search(r'Aerosol Shape \S+ Logarithmic Gaussian for Coefficient 1', svname)):
            return 0.01
        elif(re.search(r'Aerosol Shape \S+ Logarithmic Gaussian for Coefficient 2', svname)):
            return 0.001
        elif(re.search(r'Aerosol Shape \S+ Logarithmic Gaussian for Coefficient 3', svname)):
            return 1e-5
        elif(re.search(r'Zero offset waveform scale factor', svname)):
            return 0.01
        elif(re.search(r'EOF order 1 scale factor', svname)):
            return 0.01
        elif(re.search(r'Fluorescence Surface Coefficient 1', svname)):
            # This is either really small (for gosat) or really big
            # (for OCO). We use sigma_alpha to separate the 2 cases
            if(self.sigma_alpha[i] > 1e10):
                return 1e16
            else:
                return 1e-10
        elif(re.search(r'Fluorescence Surface Coefficient 2', svname)):
            return 1e-4
        # If we don't recognize name, take a stab at returning a reasonable 
        # value
        return np.sqrt(self.r.state_vector.state_covariance.diagonal()[i]) * 0.1


    @cached_property
    def apriori_covariance(self):
        '''Return the a priori covariance matrix.'''
        if(self.dist_cov_file):
            with open(self.dist_cov_file) as f:
                return pickle.load(f)
        else:
            return self.r.output_file.read_double_3d("/RetrievalResults/apriori_covariance_matrix")[0,:,:]

    @cached_property
    def measurement_uncertainty(self):
        '''Return measurement uncertainty. This is called sigma_epsilon in
        the Estimating_bias_and_Hessian.pdf document.'''
        return self.r.output_file.read_double_2d("/SpectralParameters/measured_radiance_uncert")[0,:]

    @cached_property
    def xco2_averaging_kernel(self):
        return self.r.output_file.read_double_2d("/RetrievalResults/xco2_pressure_weighting_function")[0,:]

    @cached_property
    def xco2_retrieved(self):
        return self.r.output_file.read_double_1d("/RetrievalResults/xco2")[0] * 1e6

    @property
    def xco2_portion_of_bias(self):
        '''In general, the xco2 averaging kernel only applies to the subset
        of the state vector that has the XCO2 values. This should return
        the subset of self.bias that the xco2 averaging kernel applies to.'''
        return self.bias[0:self.xco2_averaging_kernel.shape[0], 0]
