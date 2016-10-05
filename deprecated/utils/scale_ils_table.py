#!/usr/bin/env python

import os
import math
import sys
import copy
from types import ListType 

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix
from types import ListType

def scale_ils_table(input_file, output_file, scale_factor):

    # Load existing file
    print 'Reading %s' % input_file
    file_obj = OCO_Matrix(input_file)

    scale_factor = float(scale_factor)
    
    for row_idx in range(file_obj.dims[0]):
        for lbl_idx in range(file_obj.dims[1]):
            if file_obj.labels_lower[lbl_idx].find('ils_delta_lambda_') == 0:
                file_obj.data[row_idx, lbl_idx] *= scale_factor

            if file_obj.labels_lower[lbl_idx].find('ils_response_') == 0:
                file_obj.data[row_idx, lbl_idx] = (1.0/scale_factor) * file_obj.data[row_idx, lbl_idx]

    print 'Writing %s' % output_file
    file_obj.write(output_file)

def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    if type(scriptOptions) is ListType:
        scale_ils_table(fileObj.filename, fileObj.filename, scriptOptions[0])
    else:
        scale_ils_table(fileObj.filename, fileObj.filename, scriptOptions)

def standalone_main():
    if (len(sys.argv) < 4):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file> <scale_factor>\n"
        print "resample_to = <digit> or <pressure_filename>"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    scale_factor       = sys.argv[3]

    scale_ils_table(input_matrix_file, output_matrix_file, scale_factor)

if __name__ == "__main__":
    standalone_main()
