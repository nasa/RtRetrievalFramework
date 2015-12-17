import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from full_physics import *
from fluorescence import *

base_config = os.path.realpath(os.path.join(os.path.dirname(__file__), "../input/gosat/config/base_orbit_sim.lua"))
ls = LuaState.load_file(base_config)
lg = ls.globals()
lg.config = lg.BaselineConfig.new(lg.BaselineConfig)
lua_config = ls.globals().config

# Start with no flouresence
#F_COEFFS = np.array([0.0,  0.0016]) # For GOSAT data
F_COEFFS = np.array([0.0,  1.8e-3]) # For simulator data
F_FLAG = np.ones(F_COEFFS.shape[0], dtype=bool)
F_COV = np.square(np.diag([0.02, 7e-4]))

class CreatorFluorescenceEffect:
    def __init__(self, c, coeff, flag, cov):
        self.coeff = coeff
        self.flag = flag
        self.cov = cov
        
        self.c = lg.Creator.new(lg.Creator)      
        self.c.create = self.create
        self.c.initial_guess = self.initial_guess
        self.c.add_to_statevector = self.add_to_statevector

        self.f_effect = None

    def create(self, c):
        spec_index = 0

        rad = c.config.l1b.radiance(0)
        conv_units = Unit(rad.units().name()) # "Clone" to avoid lifetime issues

        self.f_effect = FluorescenceEffectPy(self.coeff, self.flag, c.config.atmosphere, c.config.l1b.sounding_zenith(0), spec_index, 13245.0, conv_units)
       
        f_output = FluorescenceOutput(self.f_effect,c.config.spec_win,c.config.instrument,conv_units)
        c.config.register_output.push_back(c.config.register_output, f_output)

        global f_plot # Have to declare global so value is available later
        f_plot = FluorescencePlot(self.f_effect)

        return (self.f_effect,)

    def initial_guess(self, c):

        # Make Fs_755 covariance as as x% of the continuum level
        cov = self.cov
        rad = c.config.l1b.radiance(0)
        cov[0,0] = np.square(np.sqrt(cov[0,0]) * np.mean(rad.data()[-10:]))
        
        igv = InitialGuessValue()
        igv.apriori_subset(self.flag, self.coeff)
        igv.apriori_covariance_subset(self.flag, cov)

        ig = CompositeInitialGuess()
        ig.add_builder(igv)

        return ig

    def add_to_statevector(self, c, sv):
        if not self.f_effect:
            raise Exception("FluorescenceEffect class not yet created")
        sv.add_observer(self.f_effect)

lua_config.fm.spectrum_effect.speceff = list(lua_config.fm.spectrum_effect.speceff) + ["fluorescence_effect"]
lua_config.fm.spectrum_effect.fluorescence_effect = { 'creator': CreatorFluorescenceEffect(lg, F_COEFFS, F_FLAG, F_COV).c }

# Write jacobians
lua_config.write_jacobian = True

# Use one window only
class CustomWindowCreator:
    def __init__(self, c):
        self.c = lg.Creator.new(lg.Creator)      
        self.c.create = self.create
        
    def create(self, c):
        win_ranges = np.array([ [(12950, 13190)],
                                [(0, 0)],
                                [(0, 0)] ])

        spec_win = SpectralWindowRange(ArrayWithUnit_double_3(win_ranges, Unit("cm^-1")))
        return spec_win

if os.environ.has_key("ABAND_ONLY"):
    lua_config.fm.spec_win.creator = CustomWindowCreator(lg).c
    lua_config.fm.atmosphere.absorber.gases = ["H2O", "O2"]

# Now create everything
lua_config.do_config(lua_config)

# Debugging stuff
if os.environ.has_key("DO_FD"):
    import finite_diff_check 
    pert_common = [ 1e-8, 1e-3 ]
    pert = [ pert_common ]
    fd_items = ["Fluorescence"]
    finite_diff_check.check_aer_coeff_fd_jac(lua_config.forward_model, fd_items, pert)

if os.environ.has_key("PLOT"):
    sv = lua_config.forward_model.state_vector()
    f_plot.solver = lua_config.solver
    sv.add_observer(f_plot)
