#!/usr/bin/env python

import os
import sys
import csv
from operator import itemgetter

import h5py
import numpy

###
# Constants
DATASET_NAME_MAPPING = {
    "Site": ("/TCCON/site_name", str),
    "sounding_id": ("/SoundingHeader/sounding_id", int),
    "xco2_true": ("/TCCON/xco2_true", float),
    "lat": ("/TCCON/latitude", float),
    "lon": ("/TCCON/longitude", float),
    "gain": ("/TCCON/gain", str),
    "surf_type": ("/TCCON/surf_type", str),
    }

SORT_BY = "/SoundingHeader/sounding_id"

VAL_FILES = sys.argv[1:]

###

tccon_data = []

for val_file in VAL_FILES:
    output_file = val_file.replace(".txt", ".h5")

    # Read input as a CSV file
    input_fobj = open(val_file, "r")

    # Extract the header
    column_headers = input_fobj.readline().strip().split()
    input_fobj.readline()

    # Read data into a list of dictionaries with the keys being the names of the HDF datasets
    # we want to create.
    for row in input_fobj.readlines():
        # Ignore rows with NaN values
        if "NaN" in row:
            continue

        data_row = {}
        for col_idx, col_val in enumerate(row.strip().split()):
            ds_name, col_type = DATASET_NAME_MAPPING[ column_headers[col_idx] ]
            data_row[ds_name] = col_type(col_val)
        tccon_data.append(data_row)

    print "%d rows read from %s" % (len(tccon_data), val_file)

    ###
    # Sort data by sounding id
    tccon_data.sort(key=itemgetter(SORT_BY))

    ###

    # Create HDF file
    out_hdf = h5py.File(output_file, 'w')

    for ds_name, ds_type in [ DATASET_NAME_MAPPING[col_name] for col_name in column_headers ]:
        g_name, ds_base = ds_name.lstrip("/").split("/")

        # Create group object if it doesnt already exist
        g_obj = out_hdf.get(g_name)
        if g_obj == None:
            g_obj = out_hdf.create_group(g_name)

        if ds_type is str:
            ds_type = h5py.special_dtype(vlen=str)

        ds_obj = g_obj.create_dataset(ds_base, (len(tccon_data),), ds_type, chunks=True)
        ds_obj.attrs['Shape'] = "Retrieval"

    # Create Datasets
    for row_idx, tccon_row in enumerate(tccon_data):
        for ds_name, ds_data in tccon_row.items():
            out_hdf[ds_name][row_idx] = ds_data

    out_hdf.close()
