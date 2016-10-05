#!/usr/bin/env python

import os
import math
import sys
import copy
from types import ListType 

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix
from types import ListType

def average_profiles(input_file_list, output_file):

    input_file_obj = open(input_file_list)
    first_file = input_file_obj.readline()
    input_file_obj.close()

    first_obj = OCO_Matrix(first_file.strip())
    dst_data = zeros((first_obj.dims[0], first_obj.dims[1]), dtype=float)
    pres_col = first_obj.labels_lower.index("pressure")
    dst_data[:, pres_col] = first_obj.data[:, pres_col]
    
    input_file_obj = open(input_file_list)

    count = 0
    for curr_atm_file in input_file_obj.readlines():
        curr_atm_file = curr_atm_file.strip()
        
        # Load existing file
        print "Loading %s" % curr_atm_file
        file_obj = OCO_Matrix(curr_atm_file)

        for col in range(file_obj.dims[1]):
            if col != pres_col:
                dst_data[:, col] += file_obj.data[:, col]

        count += 1
    
    for col in range(dst_data.shape[1]):
        if col != pres_col:        
            dst_data[:, col] /= count

    first_obj.data = dst_data
    first_obj.write(output_file)

def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    average_profiles(fileObj.filename, fileObj.filename, scriptOptions)

def standalone_main():
    if (len(sys.argv) < 3):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_atm_list_file> <output_file>\n"
        print "resample_to = <digit> or <pressure_filename>"
        sys.exit(1)

    input_atm_list_file = sys.argv[1]
    output_matrix_file  = sys.argv[2]

    average_profiles(input_atm_list_file, output_matrix_file)

if __name__ == "__main__":
    standalone_main()
