#!/usr/bin/env python

import os
import math
import sys

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

def resample_aerosol_file(input_file, output_file, output_pressure_file, verbose=False):
 
    # Load aerosol file to resample
    aero_obj = OCO_Matrix(input_file)

    # Load files containing pressure column
    pres_obj = OCO_Matrix(output_pressure_file)

    # Compute size information
    num_in_levels  = aero_obj.dims[0]
    num_out_levels = pres_obj.dims[0]

    # Load data into arrays
    press_col_in = aero_obj.labels_lower.index("pressure")
    pressures_in = numpy.zeros(num_in_levels, dtype=float)
    for row_idx in range(num_in_levels):
        pressures_in[row_idx] = aero_obj.data[row_idx][press_col_in]
        
    press_col_out = pres_obj.labels_lower.index("pressure")
    pressures_out = numpy.zeros(num_out_levels, dtype=float)
    for row_idx in range(num_out_levels):
        pressures_out[row_idx] = pres_obj.data[row_idx][press_col_out]

    # Only do resample if the grids are different since the
    # resampling operation has some issues at the top of the profile
    # and will not 100% preserve values
    if len(pressures_in) != len(pressures_out) or numpy.any(numpy.abs(pressures_in - pressures_out) > 1e-6):
        all_in_cols = range(aero_obj.dims[1])
        aero_columns = all_in_cols
        all_in_cols.remove(press_col_in)
        num_aerosols = len(aero_columns)

        if aero_obj.get_header_value('retrieval_mode').lower().find('log') == 0:
            aero_data = numpy.zeros((num_in_levels, len(aero_columns)), dtype=float)
            for col in range(len(aero_columns)):
                for row in range(num_in_levels):
                    acol = aero_columns[col]
                    aero_data[row, col] = math.exp( aero_obj.data[row, acol] )
            has_log = True
        else:
            aero_data = aero_obj.data[:, aero_columns]
            has_log = False

        ext_out = resample_aerosol(pressures_in, aero_data, pressures_out, debug=verbose)

        new_aero_data = numpy.zeros((num_out_levels, num_aerosols+1), dtype=float)
        new_aero_data[:, press_col_in] = pressures_out

        ext_idx = 0
        for aero_data_idx in aero_columns:
            if has_log:
                for row in range(num_out_levels):
                    new_aero_data[row, aero_data_idx] = math.log(ext_out[ext_idx][row])
            else:
                new_aero_data[:, aero_data_idx] = ext_out[ext_idx]
            ext_idx += 1

        aero_obj.data = new_aero_data

    # Write to output file whether we modified the data or not
    aero_obj.write(output_file)

def standalone_main():
    if (len(sys.argv) < 4):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file> <output_pressure_file>\n"
        sys.exit(1)

    input_matrix_file    = sys.argv[1]
    output_matrix_file   = sys.argv[2]
    output_pressure_file = sys.argv[3]

    resample_aerosol_file(input_matrix_file, output_matrix_file, output_pressure_file, verbose=True)

if __name__ == "__main__":
    standalone_main()
