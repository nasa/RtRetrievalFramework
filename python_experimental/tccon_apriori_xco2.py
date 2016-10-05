#!/usr/bin/env python

from full_physics import *
import sys
import os
import numpy

if len(sys.argv) < 2:
    print "Not enough args"
    sys.exit(1)

l1b_fn = sys.argv[1]
ecmwf_fn = sys.argv[2]

#print "Opening:", l1b_fn
l1b_hdf = HdfFile(l1b_fn)

sid_str = '%d' % l1b_hdf.read_double_2d("/SoundingGeometry/sounding_id")[0]
sid = OcoSoundingId(l1b_hdf, sid_str)

l1b = Level1bOco(l1b_hdf, sid)

#print "Opening:", ecmwf_fn
ecmwf = OcoSimMetEcmwf(ecmwf_fn, sid)

tccon_ap = TcconApriori(ecmwf, l1b)

static_inp = HdfFile(os.path.join(os.environ['L2_INPUT_PATH'], 'oco/input/l2_oco_static_input.h5'))
sigma_a = static_inp.read_double_1d("Pressure/Pressure_sigma_a")
sigma_b = static_inp.read_double_1d("Pressure/Pressure_sigma_b")

press = PressureSigma(sigma_a, sigma_b, ecmwf.surface_pressure, False)
press_grid = press.pressure_grid.value.value

co2_vmr = tccon_ap.co2_vmr_grid(press)
temp = TemperatureEcmwf(ecmwf, press, 0.0, False)

temp_grid = numpy.zeros(press.pressure_grid.rows, dtype=float)
for idx in range(temp_grid.shape[0]):
    p = AutoDerivativeWithUnitDouble()
    p.value = press.pressure_grid.value[idx]
    p.units = press.pressure_grid.units 
    temp_grid[idx] = temp.temperature(p).value.value

h2o_grid = ecmwf.h2o_vmr(press_grid)

from full_physics.math_util import compute_xco2
xco2 = compute_xco2(co2_vmr, press_grid, temp_grid, h2o_grid) * 1e6

print '%.2f' % xco2
