#!/usr/bin/env python

import os
import sys
import re

from OCO_Matrix import OCO_Matrix

# Quick and dirty script for outputting the control file per item jacobian files
# from the statevector jacobian file pd.dat
#
# Uses sv_names.dat from out/aggregator to figure out indexes
# and out/rad_conv.dat for start pixels header value

SV_NAMES_MATCH = { 'aero1_.*':    'aer_pd_species_1.dat',
                   'aero2_.*':    'aer_pd_species_2.dat',
                   'wc_.*':       'aer_pd_species_3.dat',
                   'ic_.*':       'aer_pd_species_4.dat',
                   'albedo_.*':   'alb_pd.dat',
                   'disp_.*':     'disp_pd.dat',
                   'co2_.*':      'mr_pd_species_1.dat',
                   'h2o_scale':   'mr_scaling_pd_species_3.dat',
                   'psurf':       'press_pd.dat',
                   'temp_offset': 't_pd_scale.dat',
                   }

def extract_sv_jacobians(pd_file, names_file, rad_conv_file):
    conv_obj = OCO_Matrix(rad_conv_file, read_data=False)
    
    names_obj = OCO_Matrix(names_file, as_strings=True)
    sv_names = names_obj['Element Name'][:,0]

    pd_obj = OCO_Matrix(pd_file)

    for name_re, output_filename in SV_NAMES_MATCH.items():
        file_indexes = []
        for curr_idx, curr_name in enumerate(sv_names):
            if re.search(name_re, curr_name):
                file_indexes.append(curr_idx)

        print output_filename, file_indexes

        out_obj = OCO_Matrix()
        out_obj.pixels = conv_obj.pixels
        out_obj.data = pd_obj.data[:,file_indexes]
        out_obj.write(output_filename)

def standalone_main():
    extract_sv_jacobians(*sys.argv[1:])

if __name__ == "__main__":
    standalone_main()
