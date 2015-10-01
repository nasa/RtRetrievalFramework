# This code was used to generate the initial table for the absorption. This
# just evaluates old SolarAbsorptionOcoFile at closely spaced points.
# This will most likely be changed in the future, but we need to be able
# to match the old calculations initially.
from full_physics import *
import numpy as np
import h5py

output_file = "solar_model_tables.h5"
if len(sys.argv) > 1:
    output_file = sys.argv[1]

print "Writing to: %s" % output_file
if os.path.exists(output_file):
    hout = h5py.File(output_file, "r+")
else:
    hout = h5py.File(output_file, "w")

s = SolarAbsorptionGfitFile(os.path.join(os.path.dirname(sys.argv[0]), "solar_merged.108"), 1.0)

band1 = np.arange(12700, 13300, 0.001)
band2 = np.arange(5800, 6500, 0.001)
band3 = np.arange(4700, 5000, 0.001)
scont1 = s.solar_absorption_spectrum(SpectralDomain(band1)).value
scont2 = s.solar_absorption_spectrum(SpectralDomain(band2)).value
scont3 = s.solar_absorption_spectrum(SpectralDomain(band3)).value

output_data = [
    ("/Solar/Absorption/Absorption_1/wavenumber", band1, "cm^-1"),
    ("/Solar/Absorption/Absorption_1/spectrum", scont1, "dimensionless"),

    ("/Solar/Absorption/Absorption_2/wavenumber", band2, "cm^-1"),
    ("/Solar/Absorption/Absorption_2/spectrum", scont2, "dimensionless"),

    ("/Solar/Absorption/Absorption_3/wavenumber", band3, "cm^-1"),
    ("/Solar/Absorption/Absorption_3/spectrum", scont3, "dimensionless"),
    ]

for ds_name, data, unit in output_data:
    if hout.get(ds_name, None):
        del hout[ds_name]
    ds = hout.create_dataset(ds_name, data=data)
    ds.attrs['Units'] = unit
    if ds_name.find("wavenumber") >= 0:
        ds.attrs['Frame'] = "Solar rest frame"

hout.close()
