#!/usr/bin/env python

import os
import math
import sys

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

def set_noise_col(input_file, output_file, snr):

    snr = float(snr)

    file_obj = OCO_Matrix(input_file)

    rad_col   = file_obj.labels_lower.index("radiance")
    noise_col = file_obj.labels_lower.index("noise")

    noise_val = max( file_obj.data[:, rad_col] ) / snr

    file_obj.data[:, noise_col] = noise_val

    file_obj.write(output_file)

def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    set_noise_col(fileObj.filename, fileObj.filename, scriptOptions)

def standalone_main():
    if (len(sys.argv) < 4):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file> <snr>\n"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    snr                = sys.argv[3]

    set_noise_col(input_matrix_file, output_matrix_file, snr)

if __name__ == "__main__":
    standalone_main()
