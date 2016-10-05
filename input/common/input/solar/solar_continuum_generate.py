# This code was used to generate the initial table for the continuum. We just
# used the old polynomial class, and made sure we had enough points that the
# old and new agreeded to better than 0.0005%. This will most likely be changed
# in the future, but we need to be able to match the old calculations initially.
from full_physics import *
import numpy as np
import scipy
from scipy.interpolate import *

def cont_calc(self, wn):
    return self.solar_continuum_spectrum(SpectralDomain([wn])).value[0]

SolarContinuumPolynomial.cont = cont_calc

hin = HdfFile("old_solar_model_input.h5", HdfFile.READ)
output_file = "solar_model_tables.h5"
if os.path.exists(output_file):
    hout = HdfFile(output_file, HdfFile.READ_WRITE)
else:
    hout = HdfFile(output_file, HdfFile.CREATE)

param = hin.read_double_with_unit_1d("/Solar/a_priori")
s = SolarContinuumPolynomial(param, False)
band1 = np.arange(12700, 13300, 0.01)
band2 = np.arange(5800, 6500, 0.01)
band3 = np.arange(4700, 5000, 0.01)
scont1 = s.solar_continuum_spectrum(SpectralDomain(band1)).value
scont2 = s.solar_continuum_spectrum(SpectralDomain(band2)).value
scont3 = s.solar_continuum_spectrum(SpectralDomain(band3)).value

ind = np.arange(12700, 13300 + 50, 50)
v = s.solar_continuum_spectrum(SpectralDomain(ind)).value
approx1 = interp1d(ind, v)
print max(abs((scont1 - approx1(band1)) / scont1 * 100))
hout.write_double_1d("/Solar/Continuum/Continuum_1/wavenumber", ind)
hout.write_double_1d("/Solar/Continuum/Continuum_1/spectrum", v)

ind = np.arange(5800, 6500 + 10, 10)
v = s.solar_continuum_spectrum(SpectralDomain(ind)).value
approx2 = interp1d(ind, v)
print max(abs((scont2 - approx2(band2)) / scont2 * 100))
hout.write_double_1d("/Solar/Continuum/Continuum_2/wavenumber", ind)
hout.write_double_1d("/Solar/Continuum/Continuum_2/spectrum", v)

ind = np.arange(4700, 5000 + 10, 10)
v = s.solar_continuum_spectrum(SpectralDomain(ind)).value
approx3 = interp1d(ind, v)
print max(abs((scont3 - approx3(band3)) / scont3 * 100))
hout.write_double_1d("/Solar/Continuum/Continuum_3/wavenumber", ind)
hout.write_double_1d("/Solar/Continuum/Continuum_3/spectrum", v)


