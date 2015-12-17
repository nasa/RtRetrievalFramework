#!/usr/bin/env python

import os
import math
import sys

from OCO_TextUtils import index_range_list
from OCO_Matrix import OCO_Matrix
from OCO_MathUtil import *

def scale_cov_by_corr(input_file, output_file, scale_factor):

    # Load existing file
    matrix_obj = OCO_Matrix(input_file)

    rows = range(matrix_obj.dims[0])
    cols = range(matrix_obj.dims[1])

    data_new = numpy.zeros((matrix_obj.dims[0], matrix_obj.dims[1]), dtype=float)

    for row_idx in rows:
        for col_idx in cols:
            rho_old = matrix_obj.data[row_idx, col_idx] / \
                      (math.sqrt(matrix_obj.data[row_idx, row_idx]) * math.sqrt(matrix_obj.data[col_idx, col_idx]))

            if rho_old < 0.0:
                sign = -1.0
            else:
                sign = 1.0
                
            fact_new = float(scale_factor) * (1.0 - abs(rho_old))
            if abs(rho_old) < 1e-40:
                rho_new = 0.0
            elif fact_new > 1.0:
                rho_new = 0.0
            else:
                rho_new = 1.0 - fact_new

            data_new[row_idx, col_idx] = sign * rho_new * \
                                         (math.sqrt(matrix_obj.data[row_idx, row_idx]) * math.sqrt(matrix_obj.data[col_idx, col_idx]))
            
    matrix_obj.data = data_new
    matrix_obj.write(output_file)
    
def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    scale_cov_by_corr(fileObj.filename, fileObj.filename, scriptOptions)


def standalone_main():
    if (len(sys.argv) < 4):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file> <columns> <offset> <method> [pressure_range]\n"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    scale_factor       = sys.argv[3]

    scale_cov_by_corr(input_matrix_file, output_matrix_file, scale_factor)

if __name__ == "__main__":
    standalone_main()
