# Short script to import JAXA update spreadsheet, which was saved as a CSV
# file
import numpy as np
import h5py

f = h5py.File("Level2/input/gosat/input/l2_gosat_static_input_quicklook.h5")
f["Instrument/ILS_jaxa_update"].create_group("ILS_1")
t = np.genfromtxt("temp.txt", skip_header=3)
a = np.empty((3,t.shape[0]))
a[0,:] = t[:,0]
a[1,:] = t[:,0]
a[2,:] = t[:,0]
g = f["Instrument/ILS_jaxa_update/ILS_1"]
g.create_dataset("delta_lambda", data=a)
# Average S and P, and also reverse order to match wavenumber
a[0,: ] = (t[:,3] + t[:,6]) / 2
a[1,: ] = (t[:,2] + t[:,5]) / 2
a[2,: ] = (t[:,1] + t[:,4]) / 2
g.create_dataset("response", data=a)
a = np.array([12900, 13050, 13200])
g.create_dataset("wavenumber", data=a)
g.create_dataset("function_type", data="interpol")
g["delta_lambda"].attrs["Units"] = "cm^-1"
g["wavenumber"].attrs["Units"] = "cm^-1"
g["response"].attrs["Units"] = "dimensionless"
