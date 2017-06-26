#!/usr/bin/env python

import os
import math
import sys

from OCO_TextUtils import index_range_list
from OCO_Matrix import OCO_Matrix
from OCO_MathUtil import *

def make_diag_only_cov(input_file, output_file):

    # Load existing file
    matrix_obj = OCO_Matrix(input_file)

    rows = range(matrix_obj.dims[0])
    cols = range(matrix_obj.dims[1])

    data_new = numpy.zeros((matrix_obj.dims[0], matrix_obj.dims[1]), dtype=float)

    for row_idx in rows:
        for col_idx in cols:
            if row_idx == col_idx:
                data_new[row_idx, col_idx] = matrix_obj.data[row_idx, col_idx]
            
    matrix_obj.data = data_new
    matrix_obj.write(output_file)
    
def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    make_diag_only_cov(fileObj.filename, fileObj.filename)


def standalone_main():
    if (len(sys.argv) < 3):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file>\n"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]

    make_diag_only_cov(input_matrix_file, output_matrix_file)

if __name__ == "__main__":
    standalone_main()
