#!/usr/bin/env python

import os
import math
import sys
import random

from OCO_TextUtils import index_range_list
from OCO_Matrix import OCO_Matrix

def random_column(input_file, output_file, columns, mean, std_dev):

    # Load existing file
    matrix_obj = OCO_Matrix(input_file)
   
    # Add ability to specify cols individually or using a * to goto end
    cols = index_range_list(columns)
    
    mean = float(mean)
    std_dev = float(std_dev)

    target_rows = range(matrix_obj.dims[0])

    for rowIdx in target_rows:
        for colIdx in cols:
            matrix_obj.data[rowIdx][colIdx] = random.gauss(mean, std_dev)

    matrix_obj.write(output_file)
    
def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    random_column(fileObj.filename, fileObj.filename, scriptOptions[0], scriptOptions[1], scriptOptions[2])


def standalone_main():
    if (len(sys.argv) < 6):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file> <columns> <mean> <std_dev>\n"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    columns            = sys.argv[3]
    mean               = sys.argv[4]
    std_dev            = sys.argv[5]

    random_column(input_matrix_file, output_matrix_file, columns, mean, std_dev)

if __name__ == "__main__":
    standalone_main()
