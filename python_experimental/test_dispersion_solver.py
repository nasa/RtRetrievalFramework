#!/usr/bin/env python
from full_physics import *
import numpy as np
import os
import sys
import math

ls = LuaState.load_file(os.path.realpath("../unit_test_data/base_config.lua"))
lg = ls.globals()
lg.config = lg.BaseConfig.new(lg.BaseConfig)
lua_config = ls.globals().config

# Instruct Lua to create everything
lua_config.do_config(lua_config)

dispersion_coefs = np.zeros((3,2), dtype=float)
dispersion_coefs[0, :] = (12869.9, 0.199493)
dispersion_coefs[1, :] = (5749.98, 0.199493)
dispersion_coefs[2, :] = (4749.93, 0.199493)

# Test old method
def solver_old_method():
    sys.path.append("../operations/populator/modules/")
    from load_disp_apriori import create_scene_dispersion_file
    from full_physics.oco_matrix import OcoMatrix

    latitude = lua_config.l1b.latitude(0).value
    sza_r = math.radians(lua_config.l1b.solar_zenith(0).value)
    saz_r = math.radians(lua_config.l1b.solar_azimuth(0).value)
    time_struct = lua_config.l1b.time().timetuple()
    aband_data = lua_config.l1b.radiance(0).data()
    apriori_out_file = "./tmp_disp_solve.dat"
    create_scene_dispersion_file(lua_config.sid_string, latitude, sza_r, saz_r, time_struct, aband_data, dispersion_coefs, apriori_out_file)
    result = np.transpose(OcoMatrix(apriori_out_file).data[:,1:])
    os.remove(apriori_out_file)
    return result

print "initial:\n", dispersion_coefs

disp_old_method = solver_old_method()
print "old method:\n", disp_old_method

print "==="
print "diff inital old:\n", disp_old_method - dispersion_coefs
print "---"

# Test that we can perturb the true solution and get back to the answer
gosat_disp = GosatDispersionFit(lua_config.l1b)
disp_new_method = gosat_disp.fit(dispersion_coefs)
print "new method sol:\n", disp_new_method
dispersion_coefs = disp_new_method.copy() 
dispersion_coefs[:, 0] += 1.0

for xx in range(3):
    gosat_disp = GosatDispersionFit(lua_config.l1b)
    disp_new_method = gosat_disp.fit(dispersion_coefs)
    print xx, "new method again:\n", disp_new_method
    print xx, "diff inital new:\n", disp_new_method - dispersion_coefs
    dispersion_coefs = disp_new_method.copy() 
