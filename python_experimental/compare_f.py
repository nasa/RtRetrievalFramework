#!/usr/bin/env python

import os
import sys
import h5py
import numpy as np
from matplotlib import pyplot as pp

l1b_file = sys.argv[1]
l2_file = sys.argv[2]

PLOT_DIR = os.path.join(os.path.dirname(l2_file), "plots", os.path.splitext(os.path.basename(l2_file))[0])

L1B_SND_DS = "/SoundingHeader/sounding_id"
L1B_WN_DS = "/SoundingHeader/wavenumber_coefficients"
L1B_F_DS = "/SoundingSpectra/fluorescence_o2"

L2_SND_DS = "/RetrievalHeader/sounding_id_reference"
L2_WN_DS = "/SpectralParameters/wavenumber_conv_fluorescence"
L2_F_DS = "/SpectralParameters/fluorescence_convolved"

if not os.path.exists(PLOT_DIR):
    os.makedirs(PLOT_DIR)

l1b_o = h5py.File(l1b_file, "r")
l2_o = h5py.File(l2_file, "r")

for idx, snd_id in enumerate(l2_o[L2_SND_DS][:]):
    l1b_idx = np.where(l1b_o[L1B_SND_DS][:] == snd_id)

    l1b_f = (np.sum(l1b_o[L1B_F_DS][l1b_idx, :, :], 1) / 2.0)[0, :]
    l1b_wn = np.poly1d(l1b_o[L1B_WN_DS][l1b_idx, 0, 0, :][0, ::-1])(range(l1b_f.shape[0]))

    l2_f = l2_o[L2_F_DS][idx, :]
    l2_wn = l2_o[L2_WN_DS][idx, :]

    pp.cla()
    pp.plot(l1b_wn, l1b_f)
    pp.plot(l2_wn, l2_f)
    pp.legend(["L1B F", "L2 F"], 0)
    pp.title("%s Fluorescence" % snd_id)
    pp.xlabel("Wavenumber (cm^-1)")

    if np.any(np.abs(l1b_f) > 0.0):
        tag = "y"
    else:
        tag = "n"

    if os.environ.has_key("SHOW"):
        pp.show()
    else:
        plot_file = os.path.join(PLOT_DIR, "f_%s_%s.png" % (snd_id, tag))
        print plot_file
        pp.savefig(plot_file)
