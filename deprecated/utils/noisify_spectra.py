#!/usr/bin/env python

import os
import sys
import math
import time
import random
import logging

from OCO_TextUtils import index_range_list 
from OCO_Matrix import OCO_Matrix
import L2_Log_Util

import numpy

random.seed(sys.argv[0])

def noisify_spectra_file(input_radiance_file, output_radiance_file, **kwarg):

    # Load existing file
    matrix_obj = OCO_Matrix(input_radiance_file)

    noisify_spectra_obj(matrix_obj, **kwargs)

    matrix_obj.write(output_radiance_file, auto_size_cols=False)    

def noisify_spectra_obj(matrix_obj, row_range_spec=None, noise_cut_off=None, save_perturb=False):
    logger = logging.getLogger(os.path.basename(__file__))

    radiance_col = matrix_obj.labels_lower.index("radiance")
    noise_col    = matrix_obj.labels_lower.index("noise")

    if row_range_spec == None or len(row_range_spec) == 0:
        row_range = range(matrix_obj.dims[0])
    else:
        row_range = index_range_list(row_range_spec)

    if matrix_obj.dims[1] <= 2 and save_perturb:
        new_data = numpy.zeros((matrix_obj.dims[0], matrix_obj.dims[1]+1), dtype=float)
        new_data[:, 0:2] = matrix_obj.data[:,:]

        matrix_obj.labels.append('True-Perturb')
        matrix_obj.units.append('W sr-1 m-2')
    else:
        new_data = matrix_obj.data

    if save_perturb:
        perturb_col  = new_data.shape[1]-1
    
    if noise_cut_off != None:
        mean_noise = 0
        for row_idx in row_range:
            mean_noise += new_data[row_idx, noise_col]
        mean_noise /= len(row_range)

    for row_idx in row_range:
        radiance_val = new_data[row_idx, radiance_col]
        noise_val    = new_data[row_idx, noise_col]

        new_radiance = random.gauss(radiance_val, noise_val)

        if (noise_cut_off != None and noise_val > float(noise_cut_off)*mean_noise):
            logger.error('Bad pixel detected for pixel row %d' % row_idx)
            noised_radiance_val = radiance_val
        else:
            noised_radiance_val = new_radiance
        
        new_data[row_idx, radiance_col] = noised_radiance_val

        if save_perturb:
            new_data[row_idx, perturb_col]  = radiance_val - noised_radiance_val
 
    matrix_obj.data = new_data

def standalone_main():
    if (len(sys.argv) < 3):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_spectra_file> <output_spectra_file> [row_range_spec]\n"
        sys.exit(1)

    input_radiance_file  = sys.argv[1]
    output_radiance_file = sys.argv[2]

    if len(sys.argv) > 3:
        row_range_spec = sys.argv[3]
    else:
        row_range_spec = None

    L2_Log_Util.init_logging()

    noisify_spectra_file(input_radiance_file, output_radiance_file, row_range_spec)

if __name__ == "__main__":
    standalone_main()
