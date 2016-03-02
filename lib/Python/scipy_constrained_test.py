from __future__ import print_function
from nose.tools import *
from full_physics import *
import scipy.optimize
from nose.plugins.skip import Skip, SkipTest

def scipy_constrained_test():
    '''This test is an example of using the scipy.optimize package. We
    don't regularly run this, and it depends on a hardcoded run. But
    leave this in place as an example of how to set this up.'''
    raise SkipTest
    config_file = "/home/smyth/Local/Level2PythonBuild/build/oco2_sounding_1_test/oco_oco2_sounding_1_test.config"
    lua_config = "/home/smyth/Local/Level2PythonBuild/build/oco2_sounding_1_test/config_diff_solv.lua"
    l2run = L2Run.create_from_existing_run(config_file, lua_config=lua_config)
    opt_problem = l2run.lua_config.opt_problem
    
    def res(x, prob):
        print("Doing res")
        print(x)
        # Note prob is smart, so if residual_x and jacobian_x is called
        # with the same x, it doesn't calculate twice
        return prob.residual_x(x)

    def res_jac(x, prob):
        print("Doing res_jac")
        print(x)
        # Note prob is smart, so if residual_x and jacobian_x is called
        # with the same x, it doesn't calculate twice
        return prob.jacobian_x(x)
    bounds = np.empty((opt_problem.parameter_size, 2))
    bounds[:,0] = -np.inf
    bounds[:,1] = np.inf
    bounds[0:20,0] = 0 # XCO2 can't go below zero
    bounds[23:(23+4*3):3, 0] = 1e-8 # Aerosol can't go below 1e-8
    # From the configuration, for GSL. These don't exactly map to 
    # scipy, but isn't far off
    max_cost_function_calls = 20
    dx_tol_rel = 1e-5
    g_tol_abs = 1e-5
    scipy.optimize.least_squares(res, opt_problem.parameters, jac=res_jac,
                                 bounds=(bounds[:,0],bounds[:,1]),
                                 max_nfev=max_cost_function_calls,
                                 xtol=dx_tol_rel,
                                 gtol=g_tol_abs,
                                 verbose=2, args=(opt_problem,))


