from config_aer_shape import *

#############

# Diagnostics for investigating aerosol retrieval behavior

# Set debug = True inside of ipython or from envrioment
from aer_diagnostics import *
sv = lua_config.forward_model.state_vector()

sv_diag = StateVectorDiagnostics(sv, lua_config.solver)

# Output SV updates in an easier to examine manner
sv.add_observer(sv_diag)

if globals().get("plot", False) or os.environ.has_key("PLOT"):
    ref_wn = 1e4/.755

    # Add observer to plot details of aerosol behavior
    plot_dir = "./plots/"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    aw = AerosolWatcher(plot_dir, lua_config.sid_string, ref_wn, lua_config.pressure, lua_config.aerosol, lua_config.solver, sv)
    sv.add_observer(aw)

# Set do_fd = True inside of ipython and %run -i <this_script>
if globals().get("do_fd", False) or os.environ.has_key("DO_FD"):
    import finite_diff_check 
    pert_common = [ 1e-8, 1e-8, 1e-8, 1e-8, 1e-8 ]
    pert = [ pert_common for x in range(4) ]
    fd_items = ["Kahn_3b Gaussian", "Kahn_2b Gaussian", "Ice Gaussian", "Water Gaussian"]
    finite_diff_check.check_aer_coeff_fd_jac(lua_config.forward_model, fd_items, pert)
