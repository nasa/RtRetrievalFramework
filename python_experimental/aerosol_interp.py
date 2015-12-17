from full_physics import *
import numpy as np

RETRIEVE = 1
HOLD_FIXED = 2
INTERPOLATE = 3

class AerosolExtinctionLogInterp(AerosolExtinctionImpBase):
    '''This is a variation of AerosolExtinctionLog that allows some of the 
    levels to be interpolated. '''
    def __init__(self, p, coeff_flag, aext, name):
        '''The flag values should be RETRIEVE, HOLD_FIXED, 
        or INTERPOLATE. aext passed in should only be for RETRIEVE and
        HOLD_FIXED values - don't pass anything in for INTERPOLATE.'''
        # Translate coeff_flag into used_flag needed by the base class.
        # The base class handles RETRIEVE and HOLD_FIXED, but not 
        # INTERPOLATE
        used_flag = []
        for f in coeff_flag:
            if(f == RETRIEVE): 
                used_flag.append(True)
            if(f == HOLD_FIXED):
                used_flag.append(False)
        self.coeff_flag = coeff_flag
        AerosolExtinctionImpBase.__init__(self,name, aext, used_flag, p)
    
    def clone(self, p):
        res = AerosolExtinctionLogInterp(p, self.coeff_flag, 
                  self.coefficient().value(), self.aerosol_name())
        return res

    def calc_aerosol_extinction(self):
        pgrid = self.pressure().pressure_grid()
        xvec = vector_auto_derivative()
        yvec = vector_auto_derivative()
        j = 0
        for i in range(pgrid.rows()):
            if(not self.coeff_flag[i] == INTERPOLATE):
                xvec.push_back(pgrid[i])
                yvec.push_back(self.coefficient()[j])
                j = j + 1
        a_inter = LinearInterpolateAutoDerivative(xvec, yvec)
        nvar = pgrid.number_variable()
        if(nvar == 0):
            nvar = self.coefficient().number_variable()
        self.aext.resize(pgrid.rows(), nvar)
        for i in range(pgrid.rows()):
            self.aext[i] = exp(a_inter(pgrid[i]))

    def state_vector_name_i(self, i):
        return "Aerosol %s Log(extinction) for Press Lvl %d" % \
            (self.aerosol_name(), i)

    def desc(self):
        coeff_flag_v = []
        for f in self.coeff_flag:
            if(f == RETRIEVE): 
                coeff_flag_v.append("RETRIEVE")
            elif(f == HOLD_FIXED):
                coeff_flag_v.append("HOLD_FIXED")
            elif(f == INTERPOLATE):
                coeff_flag_v.append("INTERPOLATE")
            else:
                coeff_flag_v.append("UNKNOWN")
        return '''
AerosolExtinction Log Interp:
   Aerosol name: %s
   Coefficient:    %s
   Coefficient flag: %s
''' % ( self.aerosol_name(), self.coefficient().value().__str__(), coeff_flag_v)

