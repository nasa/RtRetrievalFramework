#!/usr/bin/env python

import os
import sys

import h5py

base_dir = os.path.dirname(sys.argv[0])
prop_dir = os.path.join(base_dir, 'aerosol/properties/')

updated_types = sys.argv[1:]

if len(updated_types) == 0:
    print("Specify type names to update")
    sys.exit(1)

output_file = h5py.File(os.path.join(base_dir, 'l2_aerosol_combined.h5'), 'r+')

for prop_file in ['aerosol_kahn.h5', 'aerosol_merra.h5']:
    with h5py.File(os.path.join(prop_dir, prop_file), 'r') as input_file:

        for type_name in updated_types:
            if type_name in input_file.keys():
                print('%s -> %s' % (input_file.filename, type_name))

                # Delete the sub group of Properties to maintain links at the
                # top level
                src_grp = input_file[type_name]
                dst_grp = output_file[type_name]
                del dst_grp['Properties']
                output_file.copy(src_grp['Properties'], dst_grp)

output_file.close()
