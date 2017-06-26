from full_physics import *
import numpy as np

# compute total aod
def total_aod(pressure, aer_ext):
    if isinstance(pressure, np.ndarray):
        num_levels = pressure.shape[0]
    else:
        num_levels = pressure.rows()

    tot_aod = 0.0
    for layer_idx in range(num_levels-1):
        delta_press = (pressure[layer_idx + 1] - pressure[layer_idx]) / 2.0
        tot_aod +=  (delta_press * (aer_ext[layer_idx] + aer_ext[layer_idx + 1]))
    return tot_aod

class AerosolShapeGaussian(AerosolExtinctionImpBase):
    def __init__(self, pressure_obj, flag, coeffs, name, linear=False):
        ncoeff = np.array(coeffs).shape[0]
        if (ncoeff - 1) % 2 > 0:
            raise Exception("Number of coefficients must be an odd number")
        self.linear = linear
        AerosolExtinctionImpBase.__init__(self, name, coeffs, flag, pressure_obj, False)

    def clone(self, pressure_obj):
        res = AerosolShapeGaussian(pressure_obj, self.used_flag_value(), self.coefficient().value(), self.aerosol_name())
        return res

    def calc_aerosol_extinction(self, pressure_grid=None, surface_press=None):
        # Input parameters
        if self.linear:
            tot_aod = self.coefficient()[0]
        else:
            tot_aod = exp(self.coefficient()[0])

        ngaussians = int((self.coefficient().rows() - 1) / 2)

        # To allow overriding behavior in a testing setting
        if pressure_grid == None:
            pressure_grid = self.pressure().pressure_grid()
        if surface_press == None:
            surface_press = self.pressure().surface_pressure()

        if hasattr(pressure_grid, "rows"):
            number_level = pressure_grid.rows()
        else:
            number_level = len(pressure_grid)

        # We can probably add a typemap to automatically convert, but
        # for now just explictly create this
        self.aext.resize(number_level,
                         self.coefficient().number_variable())

        for g_idx in range(ngaussians):
            p0 = self.coefficient()[g_idx*2+1]
            sigma = self.coefficient()[g_idx*2+2]

            for lev in range(self.aext.rows()):
                p = pressure_grid[lev] / surface_press
                g_eval = exp( -1 * pow(p - p0, 2) / (2 * pow(sigma,2)) )

                # Because its not that easy to init ArrayAd from python to 0.0
                if g_idx == 0:
                    self.aext[lev] = g_eval
                else:
                    self.aext[lev] += g_eval
            
        curr_aod = total_aod(pressure_grid, self.aext)
        if(self.linear and tot_aod.value() < 0):
            tot_aod *= 0

        scaling_N = tot_aod/curr_aod
        for i in range(self.aext.rows()):        
            self.aext[i] = self.aext[i] * scaling_N

    def state_vector_name_i(self, i):
        return "Aerosol %s Gaussian for Coefficient %d" % (self.aerosol_name(), i)

    def desc(self):
        return '''
AerosolExtinction Gaussian:
   Coefficient:    %s
   Retrieval flag: %s
''' % ( self.coefficient().value().__str__(), self.used_flag_value().__str__())
