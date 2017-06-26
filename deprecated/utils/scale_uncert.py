#!/usr/bin/env python

import os
import sys
import math
from types import ListType

from OCO_TextUtils import index_range_list 
from OCO_Matrix import OCO_Matrix

def scale_uncertainty(input_radiance_file, output_radiance_file, scale_factor, row_range_spec=None):
       
    # Load existing file
    matrix_obj = OCO_Matrix(input_radiance_file)

    radiance_col = matrix_obj.labels_lower.index("radiance")
    noise_col    = matrix_obj.labels_lower.index("noise")

    if row_range_spec == None or len(row_range_spec) == 0:
        row_range = range(matrix_obj.dims[0])
    else:
        row_range = index_range_list(row_range_spec)

    for row_idx in row_range:
        radiance_val = matrix_obj.data[row_idx, radiance_col]
   
        new_uncert = radiance_val * float(scale_factor)
        
        matrix_obj.data[row_idx, noise_col] = new_uncert

    matrix_obj.write(output_radiance_file)

# For use by testcase generator
def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    if type(scriptOptions) is ListType:
        scale_uncertainty(fileObj.filename, fileObj.filename, scriptOptions[0], scriptOptions[1])
    else:
        scale_uncertainty(fileObj.filename, fileObj.filename, scriptOptions[0])

def standalone_main():
    if (len(sys.argv) < 4):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_spectra_file> <output_spectra_file> scale_factor [row_range_spec]\n"
        sys.exit(1)

    input_radiance_file  = sys.argv[1]
    output_radiance_file = sys.argv[2]
    scale_factor         = sys.argv[3]

    if len(sys.argv) > 4:
        row_range_spec = sys.argv[4]
    else:
        row_range_spec = None

    scale_uncertainty(input_radiance_file, output_radiance_file, scale_factor, row_range_spec)

if __name__ == "__main__":
    standalone_main()
