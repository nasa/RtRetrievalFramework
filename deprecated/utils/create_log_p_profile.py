#!/usr/bin/env python

import os
import math
import sys

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

def create_log_p_profile(input_file, output_file, column, val0, lapse_rate):

    # Load existing file
    file_obj = OCO_Matrix(input_file)
    num_rows = file_obj.dims[0]

    val0 = float(val0)
    lapse_rate = float(lapse_rate)

    # Find existing pressure bounds
    src_pres_col = file_obj.labels_lower.index("pressure")
    pressure = numpy.zeros(num_rows, dtype=float)

    for row in range(0, num_rows):
        pressure[row] = float(file_obj.data[row][src_pres_col])
   
    if column.isdigit():
        dest_prof_col = column
    else:
        dest_prof_col = file_obj.labels_lower.index(column.lower())

    # create log p profile
    for row in range(num_rows-1,0,-1):
       file_obj.data[row, dest_prof_col] = val0 - lapse_rate * (math.log(pressure[num_rows-1])-math.log(pressure[row]))

    file_obj.write(output_file)

def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    create_log_p_profile(fileObj.filename, fileObj.filename, scriptOptions[0], scriptOptions[1], scriptOptions[2])

def standalone_main():
    if (len(sys.argv) < 6):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file> <column> <val0> <lapse_rate>\n"
        print "resample_to = <digit> or <pressure_filename>"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    column             = sys.argv[3]
    val0               = sys.argv[4]
    lapse_rate         = sys.argv[5]

    create_log_p_profile(input_matrix_file, output_matrix_file, column, val0, lapse_rate)

if __name__ == "__main__":
    standalone_main()
