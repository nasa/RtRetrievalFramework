#!/usr/bin/env python
# Converts Orbit Simulator log file effective albedos
# into coefficients sufficient for L2's albedo code

import os
import sys
from optparse import OptionParser
from contextlib import nested 
import logging

import numpy
import h5py

from full_physics.acos_file import L1B
import full_physics.orbit_sim import OrbitSimLogFile

OUT_SOUNDING_ID_DS = "/SoundingGeometry/sounding_id"
OUT_ALBEDO_POLY_DS = "/Simulation/Surface/albedo_polynomial"
OUT_ALBEDO_EFF_DS = "/Simulation/Surface/eff_albedo"

logger = logging.getLogger()

def create_groups_and_dataset(hdf_obj, ds_name, *vargs, **kwargs):
    ds_parts = ds_name.strip("/").split("/")
    group_names = ds_parts[:-1]
    ds_name = ds_parts[-1]
    
    curr_group = hdf_obj
    for new_group_name in group_names:
        if new_group_name in curr_group.keys():
            curr_group = curr_group[new_group_name]
        else:
            curr_group = curr_group.create_group(new_group_name)
        
    return curr_group.create_dataset(ds_name, *vargs, **kwargs)

def extract_albedo_vals(l1b_file, log_file, output_file):

    log_obj = OrbitSimLogFile(log_file)

    with nested(L1B(l1b_file), h5py.File(output_file, "w")) as (l1b_obj, output_obj):
        sounding_ids = l1b_obj.get_sounding_ids()
        
        out_snd_id_ds = create_groups_and_dataset(output_obj, OUT_SOUNDING_ID_DS, data=sounding_ids)
        
        out_alb_poly_ds = create_groups_and_dataset(output_obj, OUT_ALBEDO_POLY_DS, dtype=float, shape=list(sounding_ids.shape) + [3,2])

        out_alb_eff_ds = create_groups_and_dataset(output_obj, OUT_ALBEDO_EFF_DS, dtype=float, shape=list(sounding_ids.shape) + [3,2])
       
        for frame_idx in range(sounding_ids.shape[0]):
            for snd_idx in range(sounding_ids.shape[1]):
                logger.info(sounding_ids[frame_idx, snd_idx])
                    
                out_alb_eff_ds[frame_idx, snd_idx, :, :] = log_obj.eff_albedo_coefficients(frame_idx, snd_idx)
                out_alb_poly_ds[frame_idx, snd_idx, :, :] = log_obj.l2_albedo_coefficients(frame_idx, snd_idx)
                

def standalone_main():
    parser = OptionParser(usage="usage: %prog <l1b_file> <log_file> -s <static_input_file>")

    parser.add_option( "-o", "--output_file", dest="output_file",
                       metavar="FILE",
                       help="destination name for extracted values")

    parser.add_option( "-v", "--verbose", dest="verbose",
                       action="store_true",
                       default=False,
                       help="enable verbose informational reporting")

    # Parse command line arguments
    (options, args) = parser.parse_args()

    # Set up logging
    if options.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    console = logging.StreamHandler()
    logger.addHandler(console)

    if len(args) < 2:
        parser.error("Need l1b and log file")

    l1b_file, log_file = args[0], args[1]

    output_file = options.output_file
    if output_file == None:
        output_file = os.path.splitext(os.path.basename(log_file))[0] + "_albedo.h5"
        
    logger.info("Reading from %s, %s output to: %s" % (l1b_file, log_file, output_file))
    extract_albedo_vals(l1b_file, log_file, output_file)

if __name__ == "__main__":
    standalone_main()
