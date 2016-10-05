from full_physics import *
import numpy as np
import math

class TemperatureEcmwf(TemperatureImpBase):
    '''This class uses ECMWF data + offset to calculate the temperature'''
    def __init__(self, ecmwf, press, temp_offset = 0.0, temp_flag = True):
        coef = np.array([temp_offset])
        flag = np.array([temp_flag])
        TemperatureImpBase.__init__(self,coef, flag, press, False)
        self.ecmwf = ecmwf

    def offset(self):
        return self.coefficient()[0]

    def offset_uncertainty(self):
        cov = self.statevector_covariance()
        if(cov.shape[0] > 0 and cov[0,0] > 0):
            return math.sqrt(cov[0,0]) 
        else:
            return 0

    def clone(self, press = None):
        if(not press):
            press = self.pressure()
        res = TemperatureEcmwf(self.ecmwf, press,
                               self.coefficient().value()[0],
                               self.used_flag_value()[0])
        return res
    
    def calc_temperature_grid(self):
        t = np.array(self.ecmwf.temperature(self.pressure().pressure_grid()))
        t += self.offset()
        self.tgrid = np_to_array_ad(t)

    def state_vector_name_i(self, i):
        return "Temperature offset (Kelvin)"

    def desc(self):
        return '''
TemperatureEcmwf:
  Temperature offset: %f
  Retrieval flag: %s
''' %(self.coefficient().value()[0], self.used_flag_value()[0].__str__()) + \
            "  Pressure: \n" + self.pressure().__str__() + \
          "  ECMWF:\n" + self.ecmwf.__str__()


class TemperatureEcmwfOutput(RegisterOutputBase):
    def __init__(self, temp):
        RegisterOutputBase.__init__(self)
        self.temp = temp

    def register_output(self, out):
        out.register_double("/RetrievalResults/temperature_offset_fph",
                            lambda : self.temp.offset().value())
        out.register_double("/RetrievalResults/temperature_offset_uncert_fph",
                            lambda : self.temp.offset_uncertainty())

    def register_output_apriori(self, out):
        frozen = self.temp.clone()
        out.register_double("/RetrievalResults/temperature_offset_apriori_fph",
                            lambda : frozen.offset().value())

# We may want to create a base class for this at some point
class CreatorTemperatureEcmwf:
    def __init__(self,lg):
        self.c = lg.Creator.new(lg.Creator)
        self.c.initial_guess = self.initial_guess
        self.c.create = self.create
        self.c.register_output = self.register_output

    def create(self,c):
        self.temp = TemperatureEcmwf(c.config.ecmwf(c.config), 
                                     c.config.pressure, 
                                     c.initial_guess(c).apriori()[0])
        return self.temp

    def initial_guess(self,c):
        read_t = c.offset
        read_t.config = c.config
        ig = InitialGuessValue()
        ig.apriori(read_t.apriori(read_t))
        ig.apriori_covariance(read_t.covariance(read_t))
        return ig
    def register_output(self,c,ro):
        ro.push_back(ro, TemperatureEcmwfOutput(self.temp))



