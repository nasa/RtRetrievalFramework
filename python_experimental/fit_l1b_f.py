#!/usr/bin/env python

import os
import sys

import h5py
import numpy as np
from config_fluorescence import *
from matplotlib import pyplot as pp
from scipy import optimize

L1B_SND_DS = "/SoundingHeader/sounding_id"
L1B_WN_DS = "/SoundingHeader/wavenumber_coefficients"
L1B_F_DS = "/SoundingSpectra/fluorescence_o2"

snd_id = 20060914120512

l1b_o = lua_config.l1b_hdf_file_v

l1b_ids = l1b_o.read_double_1d(L1B_SND_DS)
l1b_idx = np.where(l1b_ids == snd_id)

l1b_f = (np.sum(l1b_o.read_double_3d(L1B_F_DS)[l1b_idx, :, :], 2) / 2.0)[0,0,:]
l1b_wn = np.poly1d(l1b_o.read_double_4d(L1B_WN_DS)[l1b_idx, 0, 0, :][0, 0, ::-1])(range(l1b_f.shape[0]))

l1b = lua_config.l1b.radiance(0)
l1b_units = l1b.units()


def fit_func(params, plot=False):
    print >>sys.stderr, "Params:", params
    dummy_data = np.zeros(l1b_wn.shape[0])
    dummy_jac = np.zeros((l1b_wn.shape[0], len(lua_config.state_vector.state())))
    dummy_ad = ArrayAd_double_1(dummy_data, dummy_jac)
    dummy_spec = Spectrum(SpectralDomain(l1b_wn), SpectralRange(dummy_ad, l1b_units))

    f_obj = FluorescenceEffect(np.array(params), F_FLAG, lua_config.atmosphere, lua_config.l1b.sounding_zenith(0), 0, 13245.0, l1b_units)
    f_obj.apply_effect(dummy_spec)

    diff = dummy_spec.spectral_range().data() - l1b_f

    if plot:
        pp.plot(l1b_wn, l1b_f)
        pp.plot(l1b_wn, dummy_spec.spectral_range().data())
        pp.legend(["L1B", "L2"], 0)
    
    return diff

sol,success = optimize.leastsq(fit_func, list(F_COEFFS), maxfev=10000)

print "Plotting solution:", sol
fit_func(sol, plot=True)
pp.show()
