#!/usr/bin/env python

from __future__ import print_function
import sys
import os

import h5py
import numpy as np

from full_physics import *

BAND_NAME_DS = [ "o2", "weak_co2", "strong_co2" ]

# Load L2 file specified on command line
l2_fn = sys.argv[1]

if not os.path.exists(l2_fn):
    raise Exception("Could not find L2 single sounding file: %s" % l2_fn)

l2_obj = h5py.File(l2_fn, "r")

# Try and find L1B file
l1b_fn = sys.argv[2]

if not os.path.exists(l1b_fn):
    raise Exception("Could not find L1b file")

used_l1b_fn = l2_obj["/Metadata/L1BFile"][:]

# Make sure this was the one used by the L2 retrieval
#if not np.all(os.path.basename(l1b_fn) != used_l1b_fn):
#    raise Exception("Passed L1b: %s not the same as used for retrieval: %s" % (l1b_fn, used_l1b_fn))

l1b_obj = L1B(l1b_fn)

# Load L2 static input
inst_name = l1b_obj.instrument_name
static_fn = os.path.join(os.environ["L2_INPUT_PATH"], "%s/input/l2_%s_static_input.h5" % (inst_name, inst_name))

static_obj = h5py.File(static_fn, "r")
spec_win_ds = static_obj["/Spectral_Window/microwindow"]
spec_wins = ArrayWithUnit_double_3()
spec_wins.data = spec_win_ds[:]
spec_wins.units = Unit(spec_win_ds.attrs['Units'][0])

spec_win_range = SpectralWindowRange(spec_wins)
print(spec_win_range)

sounding_ids = l2_obj["/RetrievalHeader/sounding_id_reference"][:]

# Store indexes
sounding_rad_indexes = np.zeros((len(sounding_ids), len(BAND_NAME_DS), 2), dtype=int)

uniq_indexes = np.zeros((len(BAND_NAME_DS), 2), dtype=int)
# Loop over sounding ids in L2 file
for idx, sounding_id in enumerate(sounding_ids):
    shown_sounding_id = False
    dispersion = l1b_obj.get_sounding_info("dispersion", sounding_id, average='Polarization')
    disp_unit = inst_name == "gosat" and "cm^-1" or "micron"
    channel_counts = l1b_obj.get_channel_counts(sounding_id)
    for band_idx, band_name in enumerate(BAND_NAME_DS):
        band_disp = dispersion[band_idx, :]
        ret_offset = l2_obj.get("/RetrievalResults/dispersion_offset_%s" % band_name, None)
        if ret_offset != None:
            band_disp[0] = ret_offset[idx]            
        flags = np.zeros(band_disp.shape, dtype=bool)
        disp_obj = DispersionPolynomial(band_disp, flags, Unit(disp_unit), band_name, channel_counts[band_idx], True)
        indexes = spec_win_range.grid_indexes(disp_obj.pixel_grid(), band_idx)
        sounding_rad_indexes[idx, band_idx, 0] = indexes[0]
        sounding_rad_indexes[idx, band_idx, 1] = indexes[-1] 
        if(uniq_indexes[band_idx, 0] != sounding_rad_indexes[idx, band_idx, 0] or
            uniq_indexes[band_idx, 1] != sounding_rad_indexes[idx, band_idx, 1]):
            if(band_idx == 0 or not shown_sounding_id):
                print(sounding_id)
                shown_sounding_id = True
            print("  %d: (%d, %d)" % tuple([band_idx] + list(sounding_rad_indexes[idx, band_idx, :])))
            uniq_indexes[band_idx, :] = sounding_rad_indexes[idx, band_idx, :]

output_file = "l2_l1b_pixels_used.h5"
if len(sys.argv) >= 4:
    output_file = sys.argv[3]

print("Writing: %s" % output_file)
out_obj = h5py.File(output_file, "w")
out_obj.create_dataset("source_l1b_file", data=os.path.realpath(l1b_fn))
out_obj.create_dataset("source_l2_file", data=os.path.realpath(l2_fn))
out_obj.create_dataset("l1b_pixels_used", data=sounding_rad_indexes)
out_obj.create_dataset("sounding_ids", data=sounding_ids)
