from __future__ import print_function
from __future__ import division
from past.utils import old_div
from nose.tools import *
from full_physics import *
import math
try:
    import nlopt
except ImportError:
    pass                        # Ignore, these tests are just examples
from nose.plugins.skip import Skip, SkipTest

def nlopt_constrained_test():
    '''This test is an example of using the nlopt package. We ended up
    *not* using this, the solvers didn't work well with our problems. But
    leave this test in place as an example of how to set this up in case
    we want to return to this at some point.'''
    raise SkipTest
    config_file = "/home/smyth/Local/Level2PythonBuild/build/oco2_sounding_1_test/oco_oco2_sounding_1_test.config"
    lua_config = "/home/smyth/Local/Level2PythonBuild/build/oco2_sounding_1_test/config_diff_solv.lua"
    l2run = L2Run.create_from_existing_run(config_file, lua_config=lua_config)
    opt_problem = l2run.lua_config.opt_problem
    def func(x, grad):
        print(x)
        # Note gradient_x and cost_x are smart, so if called with the
        # same x it doesn't calculate twice
        if(grad.size > 0):
            g = opt_problem.gradient_x(x)
            grad[:] = g
        return opt_problem.cost_x(x)
    opt = nlopt.opt(nlopt.LD_MMA, opt_problem.parameter_size)
    opt.set_min_objective(func)
    low_x = np.zeros((opt_problem.parameter_size))
    # Not sure about these, we'll want to look at these ranges
    low_x[0:20] = 0
    low_x[20] = 0
    low_x[21] = 0
    low_x[22] = -np.inf
    low_x[23] = 1e-9
    low_x[24] = 0
    low_x[25] = 0
    low_x[26] = 1e-9
    low_x[27] = 0
    low_x[28] = 0
    low_x[29] = 1e-9
    low_x[30] = 0
    low_x[31] = 0
    low_x[32] = 1e-9
    low_x[33] = 0
    low_x[34] = 0
    low_x[35] = 0
    low_x[36] = -np.inf
    low_x[37] = 0
    low_x[38] = -np.inf
    low_x[39] = 0
    low_x[40] = -np.inf
    low_x[41] = 0
    low_x[42] = -np.inf
    low_x[43] = 0
    low_x[44] = -np.inf
    low_x[45] = 0
    low_x[46] = -np.inf
    low_x[47] = -np.inf
    low_x[48] = -np.inf
    opt.set_lower_bounds(low_x)
    x = [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
   1.00000000e-01,   1.02315449e+05,  -9.00000000e-01,   9.14938609e-01,
   0.00000000e+00,   9.50000000e-01,   1.00000000e-09,   0.00000000e+00,
   9.50000000e-01,   1.00000000e-09,   0.00000000e+00,   0.00000000e+00,
   1.00000000e-09,   0.00000000e+00,   1.00000000e+00,   0.00000000e+00,
  -9.00000000e-01,   1.13937281e+00,   9.00000000e-01,   1.01898162e+00,
  -9.00000000e-01,   1.65769130e+00,  -8.99982524e-01,   6.90711580e-01,
  -8.99963735e-01,   2.94325091e+00,   9.00046938e-01,  -2.06407296e-14,
   1.80000000e-03]
    opt.optimize(opt_problem.parameters)

def nlopt_test():
    '''This is from the tutorial'''
    raise SkipTest
    def myfunc(x, grad):
        if grad.size > 0:
            grad[0] = 0.0
            grad[1] = old_div(0.5, math.sqrt(x[1]))
        return math.sqrt(x[1])

    def myconstraint(x, grad, a, b):
        if grad.size > 0:
            grad[0] = 3 * a * (a*x[0] + b)**2
            grad[1] = -1.0
        return (a*x[0] + b)**3 - x[1]

    opt = nlopt.opt(nlopt.LD_MMA, 2)
    opt.set_lower_bounds([-float('inf'), 0])
    opt.set_min_objective(myfunc)
    opt.add_inequality_constraint(lambda x,grad: myconstraint(x,grad,2,0), 1e-8)
    opt.add_inequality_constraint(lambda x,grad: myconstraint(x,grad,-1,1), 1e-8)
    opt.set_xtol_rel(1e-4)
    x = opt.optimize([1.234, 5.678])
    minf = opt.last_optimum_value()
    print("optimum at ", x[0],x[1])
    print("minimum value = ", minf)
    print("result code = ", opt.last_optimize_result())
