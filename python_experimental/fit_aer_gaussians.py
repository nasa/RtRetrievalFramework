# Seg faults happen unless this is first, probablty
# conflicts due to library loading
from scipy import optimize

import os
import sys
import pprint

from full_physics import *
from aerosol_shape import *

from matplotlib import pyplot as pp
from matplotlib.backends.backend_pdf import PdfPages

small_set_dir = os.path.join(os.path.dirname(sys.argv[0]), "../unit_test_data/tccon_small_set")
base_config = os.path.join(small_set_dir, "config/small_set_base.lua")
ls = LuaState.load_file(base_config)

lg = ls.globals()
lg.config = lg.DynamicConfig.new(lg.DynamicConfig)
lua_config = ls.globals().config

lua_config.sid_string = "20090827005603"

lua_config.do_config(lua_config)

aerosol_names = ["Kahn_2b", "Kahn_3b", "Ice", "Water"]

# Use input pressure level grid to match against input aerosol profiles
plev = lua_config.pinp.pressure_level()

pdf_out_obj = PdfPages("aerosol_gaussian_fit.pdf")

if len(sys.argv) > 1:
    num_gauss = int(sys.argv[1])
else:
    num_gauss = 2

solution = {}
for ai, a in enumerate(aerosol_names):
    bgroup = "Aerosol/" + a + "/"
    profile_apriori = np.exp(lua_config.h.read_double_1d(bgroup + "a_priori"))

    def eval_gauss(gauss_params):
        aer_gauss_flag = np.zeros(len(gauss_params))
        ag = AerosolShapeGaussian(lua_config.pressure, aer_gauss_flag, gauss_params, a)
        # Fit to the standard input pressure grid not the local one set up by the soundign
        # we are piggy backing off of
        ag.calc_aerosol_extinction(plev, plev[-1])
        return ag.aext.value()

    def fit_func(params):
        gauss_profile = eval_gauss(params)
        diff = gauss_profile - profile_apriori
        return diff

    def jac_func(params):
        pert_common = [ 1e-10, 1e-8, 1e-8 ]
        iv = eval_gauss(params)
        jac = np.zeros((iv.shape[0], len(pert_common)))
        for i, p in enumerate(params):
            p_fw = np.copy(params)
            p_fw[i] += pert_common[i]
            c_fw = eval_gauss(p_fw)

            p_bw = np.copy(params)
            p_bw[i] -= pert_common[i]
            c_bw = eval_gauss(p_bw)

            jac[:,i] = (c_fw - c_bw) / (2 * pert_common[i])
        return jac

    sol,success = optimize.leastsq(fit_func, GAUSSIAN_COEFFS[num_gauss][a], maxfev=10000)
    solution[a] = sol

    ig = GAUSSIAN_COEFFS[num_gauss][a]
    show_ig = not np.all(np.abs(sol - ig) < 1e-5)

    f = pp.figure(ai)
    pp.cla()
    if show_ig:
        pp.plot(plev*1e-2, eval_gauss(ig))
    pp.plot(plev*1e-2, profile_apriori)
    pp.plot(plev*1e-2, eval_gauss(sol))
    pp.xlabel("Pressure (hPa)")
    pp.ylabel("Extinction (1/Pa)")
    if show_ig:
        pp.legend(["ig", "profile", "gauss"],2)
    else:
        pp.legend(["Profile", "Gaussian"],2)
    pp.title(a)

    sol_mod = [ np.exp(sol[0]) ] + [ x * plev[-1] * 1e-2 for x in sol[1:] ]
    if len(sol) > 3:
        sol_txt = r'''AOD: {0:.4f}
$\mu_1$: {1:>7.2f} hPa
$\sigma_1$: {2:>7.2f} hPa
$\mu_2$: {3:>7.2f} hPa
$\sigma_2$: {4:>7.2f} hPa'''.format(*sol_mod)
    else:
        sol_txt = r'''AOD: {0:.4f}
$\mu$: {1:>7.2f} hPa
$\sigma$: {2:>7.2f} hPa'''.format(*sol_mod)

    ax = f.axes[0]
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    ax.text(0.75, 0.25, sol_txt, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)

    pdf_out_obj.savefig(f)
print "Solution:"

pretty = pprint.PrettyPrinter(indent=4)
pretty.pprint(solution)

# Save out plots
pdf_out_obj.close()
