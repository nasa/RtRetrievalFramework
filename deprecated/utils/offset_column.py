#!/usr/bin/env python

import os
import math
import sys
import random

from OCO_TextUtils import index_range_list
from OCO_Matrix import OCO_Matrix

def offset_column(input_file, output_file, columns, offset, method, pressure_range=None):

    # Load existing file
    matrix_obj = OCO_Matrix(input_file)
   
    # Add ability to specify cols individually or using a * to goto end
    cols = index_range_list(columns)

    if offset.isdigit():
        offset = float(offset)
    else:
        offset = eval(offset)

    if pressure_range != None:
        pres_col = matrix_obj.labels_lower.index("pressure")

        pres_range_arr = pressure_range.split(',')
        pres_val_beg = float(pres_range_arr[0])
        pres_val_end = float(pres_range_arr[1])

        pres_idx_beg = 0
        pres_idx_end = matrix_obj.dims[0]

        pres_column = []
        [ pres_column.append(float(val[pres_col])) for val in matrix_obj.data ]

        pres_idx_curr = 0
        beg_found = False
        for pres_val in pres_column:            
            if pres_val >= pres_val_beg and not beg_found:
                pres_idx_beg = pres_idx_curr
                beg_found = True
        
            if pres_val <= pres_val_end:
                pres_idx_end = pres_idx_curr + 1

            pres_idx_curr += 1

        target_rows = range(pres_idx_beg, pres_idx_end)

    else:
        target_rows = range(matrix_obj.dims[0])

    for rowIdx in target_rows:
        for colIdx in cols:

            #print 'old_val[%d][%d] = %f' % (rowIdx, colIdx, matrix_obj.data[rowIdx][colIdx])
            
            if method == '/':
                matrix_obj.data[rowIdx][colIdx] = matrix_obj.data[rowIdx][colIdx] / offset
            elif method == '-':
                matrix_obj.data[rowIdx][colIdx] = matrix_obj.data[rowIdx][colIdx] - offset
            elif method == '*':
                matrix_obj.data[rowIdx][colIdx] = matrix_obj.data[rowIdx][colIdx] * offset
            else:
                matrix_obj.data[rowIdx][colIdx] = matrix_obj.data[rowIdx][colIdx] + offset

            #print 'new_val[%d][%d] = %f' % (rowIdx, colIdx, matrix_obj.data[rowIdx][colIdx])

    matrix_obj.write(output_file)
    
def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    if len(scriptOptions) > 3:
        offset_column(fileObj.filename, fileObj.filename, scriptOptions[0], scriptOptions[1], scriptOptions[2], scriptOptions[3])
    else:
        offset_column(fileObj.filename, fileObj.filename, scriptOptions[0], scriptOptions[1], scriptOptions[2])


def standalone_main():
    if (len(sys.argv) < 6):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file> <columns> <offset> <method> [pressure_range]\n"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    columns            = sys.argv[3]
    offset             = sys.argv[4]
    method             = sys.argv[5]

    if len(sys.argv) > 6:
        pressure_range = sys.argv[6]
    else:
        pressure_range = None

    offset_column(input_matrix_file, output_matrix_file, columns, offset, method, pressure_range)

if __name__ == "__main__":
    standalone_main()
