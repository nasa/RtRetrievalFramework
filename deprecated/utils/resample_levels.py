#!/usr/bin/env python

import os
import math
import sys
import copy
from types import ListType 

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix
from types import ListType

def resample_levels(input_file, output_file, resample_to, val_extrapolate=False):

    # Load existing file
    file_obj = OCO_Matrix(input_file)

    try:
        src_pres_col = file_obj.labels_lower.index("pressure")
    except:
        raise IOError('Could not find pressure column in input file: "%s"' % input_file)

    try:
        src_temp_col = file_obj.labels_lower.index("t")
    except:
        src_temp_col = -1

    ## Do nothing except write output file if input and desired levels already match
    if resample_to.isdigit() and file_obj.dims[0] == int(resample_to):
        file_obj.write(output_file)
        return

    elif resample_to.isdigit():
        resample_to = int(resample_to)
        dst_data = numpy.zeros((resample_to, file_obj.dims[1]), dtype=float)
        
    elif os.path.exists(resample_to):
        dest_pressure_file = resample_to
        pres_obj = OCO_Matrix(dest_pressure_file)

        dst_pres_col = pres_obj.labels_lower.index("pressure")

        dst_data = numpy.zeros((pres_obj.dims[0], file_obj.dims[1]), dtype=float)

        resample_to = pres_obj.data[:, dst_pres_col]
        
    else:
        raise ValueError('Resample to argument "%s" is neither an integer nor a file that exists' % resample_to)


    for col_idx in range(file_obj.dims[1]):
        # Interpolate all but temperature in log space
        if col_idx == src_temp_col:
            log_data=False
        else:
            log_data=True

        if col_idx == src_pres_col:
            do_extrapolate = True
        else:
            do_extrapolate = val_extrapolate
            
        dst_data[:, col_idx] = resample_profile( file_obj.data[:, src_pres_col],
                                                 file_obj.data[:, col_idx], 
                                                 resample_to,
                                                 log_data=log_data,
                                                 extrapolate=do_extrapolate )
    file_obj.data = dst_data
    file_obj.write(output_file)

def standalone_main():
    if (len(sys.argv) < 4):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file> <resample_to> [<extrapolate_t_f>]\n"
        print "resample_to = <digit> or <pressure_filename>"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    resample_to        = sys.argv[3]

    if len(sys.argv) > 4:
        extrapolate    = bool(sys.argv[4])
    else:
        extrapolate    = False

    resample_levels(input_matrix_file, output_matrix_file, resample_to, extrapolate)

if __name__ == "__main__":
    standalone_main()
