from full_physics import *
import numpy as np
import math

RETRIEVE = 1
HOLD_FIXED = 2
INTERPOLATE = 3

class AbsorberVmrLevel(AbsorberVmrImpBase):
    '''This class just returns the cofficients as the vmr at each pressure level. '''
    def __init__(self, gas_name, vmr, coeff_flag, press):
        '''Initialize. The flag values should be RETRIEVE, HOLD_FIXED, 
        or INTERPOLATE. vmr passed in should only be for RETRIEVE and
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
        AbsorberVmrImpBase.__init__(self, gas_name, vmr, used_flag, press, 
                                    False)

    def clone(self, press = None):
        if(not press):
            press = self.pressure()
        res = AbsorberVmrLevel(self.gas_name(), 
                               self.coefficient().value(),
                               self.coeff_flag,
                               press)
        return res
    
    def calc_vmr(self):
        pgrid = self.pressure().pressure_grid()
        xvec = vector_auto_derivative()
        yvec = vector_auto_derivative()
        j = 0
        for i in range(pgrid.rows()):
            if(not self.coeff_flag[i] == INTERPOLATE):
                xvec.push_back(pgrid[i])
                yvec.push_back(self.coefficient()[j])
                j = j + 1
        vmr_inter = LinearInterpolateAutoDerivative(xvec, yvec)
        nvar = pgrid.number_variable()
        if(nvar == 0):
            nvar = self.coefficient().number_variable()
        self.vmr.resize(pgrid.rows(), nvar)
        for i in range(pgrid.rows()):
            self.vmr[i] = vmr_inter(pgrid[i])

    def state_vector_name_i(self, i):
        return "%s VMR for Press Lvl %d" %(self.gas_name(), i + 1)

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
AbsorberVmrLevel:
  Gas: %s
  VMR Coefficient: %s
  Coefficient flag: %s
''' %(self.gas_name(), self.coefficient().value().__str__(), 
      coeff_flag_v)

class AbsorberVmrEcmwf(AbsorberVmrImpBase):
    '''This class just returns the cofficients as the vmr at each pressure level'''
    def __init__(self, gas_name, vmr, flag, press, ecmwf):
        AbsorberVmrImpBase.__init__(self, gas_name, vmr, flag, press, False)
        self.ecmwf = ecmwf

    def scale(self):
        return self.coefficient()[0]

    def scale_uncertainty(self):
        cov = self.statevector_covariance()
        if(cov.shape[0] > 0 and cov[0,0] > 0):
            return math.sqrt(cov[0,0]) 
        else:
            return 0

    def clone(self, press = None):
        if(not press):
            press = self.pressure()
        res = AbsorberVmrEcmwf(self.gas_name(), 
                               self.coefficient().value(),
                               self.used_flag_value(),
                               press, self.ecmwf)
        return res
    
    def calc_vmr(self):
        t = np.array(self.ecmwf.h2o_vmr(self.pressure().pressure_grid())) * \
            self.scale()
        self.vmr = np_to_array_ad(t)

    def state_vector_name_i(self, i):
        return "%s Scaling factor" % self.gas_name()

    def desc(self):
        return '''
AbsorberVmrEcmwf:
  Gas:            %s
  Scale:          %s
  Retrieval flag: %s
  ECMWF:          %s
''' %(self.gas_name(), self.scale().value().__str__(), 
      self.used_flag_value()[0].__str__(), self.ecmwf.__str__()) 

class AbsorberVmrEcmwfOutput(RegisterOutputBase):
    def __init__(self, abs):
        RegisterOutputBase.__init__(self)
        self.abs = abs

    def register_output(self, out):
        out.register_double("/RetrievalResults/h2o_scale_factor",
                            lambda : self.abs.scale().value())
        out.register_double("/RetrievalResults/h2o_scale_factor_uncert",
                            lambda : self.abs.scale_uncertainty())

    def register_output_apriori(self, out):
        frozen = self.abs.clone()
        out.register_double("/RetrievalResults/h2o_scale_factor_apriori",
                            lambda : frozen.scale().value())


class CreatorAbsorberVmrLevel:
    def __init__(self, lg, hold_fixed = False):
        self.c = lg.CreatorVmr.new(lg.CreatorVmr)
        self.c.create_vmr = self.create_vmr
        self.c.initial_guess = self.initial_guess
        self.hold_fixed = hold_fixed

    def initial_guess(self, c):
        if(self.hold_fixed):
            return CompositeInitialGuess()
        res = InitialGuessValue()
        res.apriori(c.apriori_v(c))
        res.apriori_covariance(c.covariance_v(c))
        return res

    def create_vmr(self, c):
        a = c.apriori_v(c)
        coeff = np.empty(a.size)
        if(self.hold_fixed):
            coeff[:] = HOLD_FIXED
        else:
            coeff[:] = RETRIEVE
        return AbsorberVmrLevel(c.name, c.apriori_v(c), coeff, 
                                c.config.pressure)

class CreatorAbsorberVmrEcmwf:
    def __init__(self, lg):
        self.c = lg.ConfigCommon.vmr_fixed_level_scaled.new(lg.ConfigCommon.vmr_fixed_level_scaled)
        self.c.create_vmr = self.create_vmr

    def create_vmr(self, c):
        return AbsorberVmrEcmwf(c.name, [c.scale_apriori],
                               [True], c.config.pressure, 
                                c.config.ecmwf(c.config))

def create_absorber_ecmwf(c):
    vmr = vector_absorber_vmr()
    absco = vector_gas_absorption()
    cig = CompositeInitialGuess()
    read_gas = c.creators.gas
    # CO2
    apriori = read_gas.CO2.apriori(c)
    absorb = AbsorberVmrLevel("CO2", apriori, read_gas.CO2.coeff_flag, 
                              c.pressure)
    vmr.push_back(absorb)
    ig = InitialGuessValue()
    apriori_flag = absorb.used_flag_value()
    ig.apriori_subset(apriori_flag, apriori)
    ig.apriori_covariance_subset(apriori_flag, read_gas.CO2.covariance(c))
    cig.add_builder(ig)
    absco.push_back(c.open_absco(c, read_gas.CO2.absco, 1.0))

    #H2O
    h2o_scale_apriori = [read_gas.H2O.scale_apriori]
    h2o_scale_cov = [[read_gas.H2O.scale_cov]]
    h2o_vmr = AbsorberVmrEcmwf("H2O", h2o_scale_apriori,
                               [True], c.pressure, c.ecmwf)
    vmr.push_back(h2o_vmr)
    ig = InitialGuessValue()
    ig.apriori(h2o_scale_apriori)
    ig.apriori_covariance(h2o_scale_cov)
    cig.add_builder(ig)
    absco.push_back(c.open_absco(c, read_gas.H2O.absco, 1.0))

    #O2
    apriori = read_gas.O2.apriori(c)
    coeff_flag = []
    for i in range(apriori.size):
        coeff_flag.append(HOLD_FIXED)
    vmr.push_back(AbsorberVmrLevel("O2", apriori, coeff_flag, c.pressure))
    absco.push_back(c.open_absco(c, read_gas.O2.absco, 1.0))

    c.absorber = AbsorberAbsco(vmr, c.pressure, c.temperature,
                               c.altitude, absco)
    c.ab_ig = cig

    # Output
    c.register_output.push_back(c.register_output,
                                AbsorberVmrEcmwfOutput(h2o_vmr))
    
