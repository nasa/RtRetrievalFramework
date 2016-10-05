from full_physics import *
import numpy as np

class PressureSigma(PressureImpBase):
    '''This class uses pressure sigma levels to determine the pressure
    levels we do the RT on.'''
    def __init__(self, a, b, surface_pressure, pressure_flag = True):
        coef = np.array([surface_pressure])
        flag = np.array([pressure_flag])
        PressureImpBase.__init__(self,coef, flag)
        self.a = a
        self.b = b

    def clone(self):
        res = PressureSigma(self.a, self.b,
                            self.coefficient().value()[0],
                            self.used_flag_value()[0])
        return res
    
    def calc_pressure_grid(self):
        t = self.b * self.coefficient()[0] + self.a
        self.pgrid.reference(np_to_array_ad(t))

    def state_vector_name_i(self, i):
        return "Surface Pressure (Pascals)"

    def desc(self):
        return '''
Pressure Sigma:
  Surface pressure: %f
  Retrieval flag: %s
  a: %s
  b: %s
''' %(self.coefficient().value()[0], self.used_flag_value()[0].__str__(),
      self.a.__str__(), self.b.__str__())

class CreatorPressureSigma:
    def __init__(self, lg, a, b):
        self.c = lg.CreatorApriori.new(lg.CreatorApriori)
        self.c.create = self.create
        self.a = a
        self.b = b

    def create(self,c):
        return PressureSigma(self.a, self.b, c.initial_guess(c).apriori()[0])
        

