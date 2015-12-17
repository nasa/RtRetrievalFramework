from full_physics import *
import numpy as np

import os
import sys
from functools import wraps
from matplotlib import pyplot

def get_aer_coeff_idx(sv, aer_names):
    aer_sv_idx = []
    for nm in aer_names:
        coeff_idx = [y[0] for y in filter(lambda x: x[1].find("%s Polynomial" % nm) >= 0, enumerate(sv.state_vector_name()))]
        aer_sv_idx.append(coeff_idx)
    return aer_sv_idx

class StateVectorDiagnostics(ObserverStateVector):
    def __init__(self, state_vector, solver):
        ObserverStateVector.__init__(self)
        self.state_vector = state_vector
        self.solver = solver
        self.prev_sv = np.copy(state_vector.state())
        
    def notify_add(self, sv):
        pass

    def notify_remove(self, sv):
        pass

    def notify_update(self, sv):
        sv = self.state_vector
        sv_name = sv.state_vector_name()
        sv_curr = sv.state()

        for prev, curr, name in zip(self.prev_sv, sv_curr, sv_name):
            print "{0:>40} {1:>18} {2:>18}".format(name[:40], curr, (curr-prev))
        self.prev_sv = np.copy(sv_curr)

        iteration = self.solver.number_iteration()
        divergence = self.solver.number_divergent()
        print "\nIteration #%d, Divergence #%d\n" % (iteration, divergence)

#class AerosolWatcher(ObserverConnorSolver):
class AerosolWatcher(ObserverStateVector):
    def __init__(self, plot_dir, sounding_id, ref_wn, pressure, aerosol, solver, state_vector):
        ObserverStateVector.__init__(self)
        self.plot_dir = plot_dir
        self.sounding_id = sounding_id
        self.ref_wn = ref_wn
        self.pressure = pressure
        self.aerosol = aerosol
        self.solver = solver
        self.state_vector = state_vector
        
    def notify_add(self, solver):
        pass
    def notify_remove(self, solver):
        pass

    def plot_aer_od(self, iteration, divergence):
        aer_names = self.aerosol.aerosol_name()
        aer_od_val = self.aerosol.optical_depth_each_layer(self.ref_wn).value()
        aer_od_jac = self.aerosol.optical_depth_each_layer(self.ref_wn).jacobian()

        pyplot.figure(0)
        pyplot.cla()

        press_levels = self.pressure.pressure_grid().value()
        press_layers = (press_levels[:-1] + press_levels[1:]) / 2.0
        press_diff = (press_levels[1:] - press_levels[:-1])

        for aer_idx in range(aer_od_val.shape[1]):
            pyplot.plot(press_layers*1e-2, aer_od_val[:,aer_idx])
        pyplot.legend(aer_names,0)
        pyplot.title("Iteration %d Divergence %d for Wn: %.2f" % (iteration, divergence, self.ref_wn))
        plot_file = os.path.join(self.plot_dir, "aer_od_%s_i%02d_d%02d.png" % (self.sounding_id, iteration, divergence))
        pyplot.savefig(plot_file)

    def plot_aer_coeffs(self, iteration, divergence):
        sv = self.state_vector
        aer_names = self.aerosol.aerosol_name()
        press_levels = self.pressure.pressure_grid().value()
        
        aer_coeffs = get_aer_coeff_idx(sv, aer_names)
        if len(aer_coeffs[0]) == 0:
            return

        pyplot.cla()
        coeff_aer_names = []
        for aer_idx, a_name in zip(aer_coeffs, aer_names):
            sv_val = sv.state_with_derivative().value()[aer_idx]
            if len(sv_val) > 0:
                print >>sys.stderr, "%s aerosol polynomial coeffs: %s" % (a_name, repr(sv_val))
                pyplot.plot( press_levels*1e-2, np.exp(np.lib.polynomial.poly1d(sv_val)(press_levels)) )
                coeff_aer_names.append(a_name)
        pyplot.legend(coeff_aer_names,0)
        pyplot.title("Aer Coeffs Iteration %d Divergence %d\n" % (iteration, divergence))
        plot_file = os.path.join(self.plot_dir, "aer_coeff_%s_i%02d_d%02d.png" % (self.sounding_id, iteration, divergence))
        pyplot.savefig(plot_file)

    def notify_update(self, solver):
        iteration = self.solver.number_iteration()
        divergence = self.solver.number_divergent()
        print >>sys.stderr, "Solver iteration: %d, div: %d, plotting aerosol stuff" % (iteration, divergence)
        self.plot_aer_od(iteration, divergence)
