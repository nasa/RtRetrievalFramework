# This code was used to generate the initial table for the absorption. This
# just evaluates old SolarAbsorptionOcoFile at closely spaced points.
# This will most likely be changed in the future, but we need to be able
# to match the old calculations initially.
from full_physics import *
import numpy as np


h = HdfFile("/home/drt/acos/development/input/fts/input/l2_fts_static_input.h5", HdfFile.READ_WRITE)

frac = 0.0024 / 0.0092 # the field of view FOVO parameter, divided by average solar disk size

s = SolarAbsorptionOcoFile(h, "Solar", frac)
band1 = np.arange(12700, 13300, 0.01)
band2 = np.arange(5800, 6500, 0.01)
band3 = np.arange(4700, 5000, 0.01)
scont1 = s.solar_absorption_spectrum(SpectralDomain(band1)).value
scont2 = s.solar_absorption_spectrum(SpectralDomain(band2)).value
scont3 = s.solar_absorption_spectrum(SpectralDomain(band3)).value

h.write_double_1d("/Solar/Absorption/Absorption_1/wavenumber", band1)
h.write_double_1d("/Solar/Absorption/Absorption_1/spectrum", scont1)

h.write_double_1d("/Solar/Absorption/Absorption_2/wavenumber", band2)
h.write_double_1d("/Solar/Absorption/Absorption_2/spectrum", scont2)

h.write_double_1d("/Solar/Absorption/Absorption_3/wavenumber", band3)
h.write_double_1d("/Solar/Absorption/Absorption_3/spectrum", scont3)


