#!/usr/bin/env python

import os
import math
import sys
import copy
import logging

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix
from types import ListType
from OCO_TextUtils import index_range_list
import L2_Log_Util

def create_simple_cov(input_file, output_file, scaling, columns):
    logger = logging.getLogger(os.path.basename(__file__))

    # Load source file
    src_obj = OCO_Matrix(input_file)

    scaling = [float(v) for v in scaling.split(',')]
    if len(scaling) < src_obj.dims[0]:
        last_val = scaling[len(scaling)-1]
        [ scaling.append(last_val) for nada in range(src_obj.dims[0] - len(scaling)) ]

    try:
        columns = index_range_list(columns, max_value=src_obj.dims[1])
    except:
        if not type(columns) is ListType:
            col_name_list = [columns]
        else:
            col_name_list = columns

        columns = []
        for curr_name in col_name_list:
            if curr_name.lower() not in src_obj.labels_lower:
                raise IOError('Column named %s not found in file: %s' % (curr_name, input_file))
            columns.append( src_obj.labels_lower.index(curr_name.lower()) )
                                    
    logger.info('cols = ', columns)

    num_diags = len(columns) * src_obj.dims[0]
    logger.info('num_diags = ', num_diags)
    dst_data = numpy.zeros((num_diags, num_diags), dtype=float)

    diag_index = 0
    for col_index in columns:
        for row_index in range(src_obj.dims[0]):
            #print '%d, %d => %d, %d' % (row_index, col_index, diag_index, diag_index)
            dst_data[diag_index, diag_index] = (src_obj.data[row_index, col_index] * scaling[row_index])**2
            diag_index += 1

    logger.info('Writing: %s' % input_file)
    src_obj.file_id = 'Simple Covariance Created from "%s", scaling: "%s"' % (input_file, ', '.join([str(sc) for sc in scaling]))
    src_obj.labels = []
    src_obj.data = dst_data
    src_obj.units = []
    src_obj.write(output_file, auto_size_cols=False, verbose=True)

def standalone_main():
    if (len(sys.argv) < 5):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_covariance_file> <scaling> <column_def>\n"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    scaling            = sys.argv[3]
    column_def         = sys.argv[4]

    L2_Log_Util.init_logging()
    create_simple_cov(input_matrix_file, output_matrix_file, scaling, column_def)

if __name__ == "__main__":
    standalone_main()
