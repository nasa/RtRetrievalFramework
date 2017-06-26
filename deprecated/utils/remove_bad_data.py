#!/usr/bin/env python

import re
import os
import math
import sys

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

def remove_bad_data_all(input_file, output_file, check_col, check_val):

    # Load existing file
    file_obj = OCO_Matrix(input_file)
    num_rows = file_obj.dims[0]

    if check_col.isdigit():
        check_col = int(check_col)
    else:
        check_col = file_obj.labels_lower.index(check_col.lower())
    
    good_mask = []

    for row_idx in range(num_rows):
        if not re.search(str(check_val).lower(), str(file_obj.data[row_idx, check_col]).lower()):
            good_mask.append(row_idx)

    cleaned_data = numpy.zeros((len(good_mask), file_obj.dims[1]), dtype=float)

    new_data_idx = 0
    for good_row in good_mask:
        cleaned_data[new_data_idx, :] = file_obj.data[good_row, :]
        new_data_idx += 1
    
    file_obj.data = cleaned_data
    file_obj.write(output_file)

def remove_bad_data_last(input_file, output_file, check_col, check_val):

    # Load existing file
    file_obj = OCO_Matrix(input_file)
    num_rows = file_obj.dims[0]

    if check_col.isdigit():
        check_col = int(check_col)
    else:
        check_col = file_obj.labels_lower.index(check_col.lower())

    last_good_index = -1
    for row_idx in range(num_rows-1, 1, -1):
        if not re.search(str(check_val).lower(), str(file_obj.data[row_idx, check_col]).lower()):
            last_good_index = row_idx
            break
            
    print "Last good index = ", last_good_index

    file_obj.dims = [last_good_index+1, file_obj.dims[1]]
    file_obj.write(output_file, use_set_dims=True, auto_size_cols=False)

def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    remove_bad_data_last(fileObj.filename, fileObj.filename, scriptOptions[0], scriptOptions[1])

def standalone_main():
    if (len(sys.argv) < 5):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_matrix_file> <check_col> <check_val>\n"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    check_col          = sys.argv[3]
    check_val          = sys.argv[4]

    remove_bad_data_last(input_matrix_file, output_matrix_file, check_col, check_val)

if __name__ == "__main__":
    standalone_main()
