from full_physics import *
import numpy as np

class AerosolExtinctionLog10(AerosolExtinctionImpBase):
    '''This is a class that uses Log10 instead of Log distribution. This is
    just a toy example, there isn't any advantage of this over the
    AerosolExtinctionLog we already have. But this is a good test of all
    the features of writing in python'''
    def __init__(self, p, flag, aext, name):
        AerosolExtinctionImpBase.__init__(self,name, aext, flag, p)
    
    def clone(self, p):
        res = AerosolExtinctionLog10(p, self.used_flag_value(), 
                  self.coefficient().value(), self.aerosol_name())
        return res

    def calc_aerosol_extinction(self):
        coef = np.array(self.coefficient())
        self.aext = np_to_array_ad(pow(10, coef))

    def state_vector_name_i(self, i):
        return "Aerosol %s Log10(extinction) for Press Lvl %d" % \
            (self.aerosol_name(), i)

    def desc(self):
        return '''
AerosolExtinction Log 10:
   Coefficient:    %s
   Retrieval flag: %s
''' % ( self.coefficient().value().__str__(), self.used_flag_value().__str__())

def create_aerosol_log_10(c):
    aer_flag = np.ones(c.number_pressure_level, np.bool)
    vap = vector_aerosol_property();
    vex = vector_aerosol_extinction();
    # We can try reading this from lua a bit later, but right now we don't
    # directly support reading arrays.
    aerosols = ["Kahn_2b", "Kahn_3b", "Water", "Ice"]
    cig = CompositeInitialGuess()
    c.number_aerosol = len(aerosols)
    for a in aerosols:
        bgroup = "Aerosol/" + a + "/"
        apriori = c.h.read_double_1d(bgroup + "a_priori")
        cov = c.h.read_double_2d(bgroup + "covariance")
        # convert to log10
        apriori = np.log10(np.exp(apriori))
        cov = np.square(np.log10(np.exp(np.sqrt(cov))))
        vap.push_back(AerosolPropertyHdf(c.h, bgroup + "Properties"))
        vex.push_back(AerosolExtinctionLog10(c.pressure, aer_flag, 
                                             apriori, a))
        ig = InitialGuessValue()
        ig.apriori_subset(aer_flag, apriori)
        ig.apriori_covariance_subset(aer_flag, cov)
        cig.add_builder(ig)
    c.aerosol = Aerosol(vex, vap, c.pressure)
    c.a_ig = cig

