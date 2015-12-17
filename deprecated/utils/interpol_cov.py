#!/usr/bin/env python

import os
import math
import sys

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

def interpol_cov(input_file, output_file, src_pressure_file, dst_pressure_file):

    # Load existing file
    file_obj = OCO_Matrix(input_file)

    # Find existing pressure bounds
    pres_src_obj = OCO_Matrix(src_pressure_file)
    src_pres_col = pres_src_obj.labels_lower.index("pressure")
    num_levels_src = pres_src_obj.dims[0]
    pressure_src = pres_src_obj.data[:, src_pres_col]

    pres_dst_obj = OCO_Matrix(dst_pressure_file)
    dst_pres_col = pres_dst_obj.labels_lower.index("pressure")
    num_levels_dst = pres_dst_obj.dims[0]
    pressure_dst = pres_dst_obj.data[:, dst_pres_col]

    M = numpy.zeros((num_levels_dst, num_levels_src), dtype=float)

    # Setup Interpolation Matrix 
    for i in range(num_levels_dst):
        for j in range(num_levels_src):
            if (pressure_dst[i] <= pressure_src[j]):
                lev = j
                break

        if (lev > 0):
            M[i, lev] = \
                 (math.log(pressure_dst[i]) - math.log(pressure_src[lev-1]))      \
                 / (math.log(pressure_src[lev])                                \
                    - math.log(pressure_src[lev-1]))
        else:
            M[i, lev] = 1.0

        if (i > 0):
            M[i, lev-1] =                                 \
                 (-math.log(pressure_dst[i]) + math.log(pressure_src[lev-1]))     \
                 / (math.log(pressure_src[lev])                                \
                    - math.log(pressure_src[lev-1])) + 1

    # Use interpolation matrix to create new covariance
    shat_out = (mat(M) * mat(file_obj.data)) * transpose(mat(M))

    file_obj.data = shat_out
    file_obj.write(output_file, auto_size_cols=False)

def standalone_main():
    if (len(sys.argv) < 5):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file> <src_pressure_file> <dst_pressure_file>\n"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    src_pressure_file  = sys.argv[3]
    dst_pressure_file  = sys.argv[4]

    interpol_cov(input_matrix_file, output_matrix_file, src_pressure_file, dst_pressure_file)

if __name__ == "__main__":
    standalone_main()
