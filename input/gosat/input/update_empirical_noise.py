# Short script to update emperical noise, easier than using HDFview to do this
import h5py
import numpy as np

f = h5py.File("l2_gosat_static_input.h5")
t = f["/Instrument/Empirical_Noise"][:,:]
print t
f["/Instrument/Empirical_Noise"][:,:] = np.array(
    [[0.0027280862,0.0020104083,0.0028957409],
     [0.87458300,  0.82221914, 0.80018426]])

