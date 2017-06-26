#!/usr/bin/env python

import os
import math
import sys
import logging

from Generator_Utils import *
from OCO_MathUtil import *
from OCO_TextUtils import index_range_list, evaluate_bool_str
from OCO_Matrix import OCO_Matrix

import numpy
import random

# For consistency if random values are used in 'modify' statement
random.seed(sys.argv[0])

import copy as copy_module

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    logger = logging.getLogger(os.path.basename(__file__))

    # Load existing file
    matrix_obj = OCO_Matrix(source)

    for modifySect in moduleSections:

        # Add ability to specify cols individually or using a * to goto end
        columns = Apply_Template(modifySect.Get_Keyword_Value('columns'), valuesDict, mapDict=mapDict)
        rows    = Apply_Template(modifySect.Get_Keyword_Value('rows'), valuesDict, mapDict=mapDict)
        modify  = modifySect.Get_Keyword_Value('modify')
        delete  = evaluate_bool_str( modifySect.Get_Keyword_Value('delete') )
        add_column = evaluate_bool_str( modifySect.Get_Keyword_Value('add_column') )

        if columns != None:
            try:
                columns = index_range_list(columns, max_value=matrix_obj.dims[1])
            except:
                if not type(columns) is ListType:
                    col_name_list = [columns]
                else:
                    col_name_list = columns

                columns = []
                for curr_name in col_name_list:
                    if curr_name.lower() not in matrix_obj.labels_lower:
                        if add_column:
                            matrix_obj.add_column(curr_name)
                            columns.append(  matrix_obj.dims[1] - 1 )
                        else:
                            raise IOError('Column named %s not found in file: %s' % (curr_name, source))
                            
                    columns.append( matrix_obj.labels_lower.index(curr_name.lower()) )
        else:
            columns = range(matrix_obj.dims[1])

        if rows != None:
            rows = index_range_list(rows, max_value=matrix_obj.dims[0])
        else:
            rows = range(matrix_obj.dims[0])

        if delete and modify != None:
            raise ValueError('delete and modify keywords can not be specified together')

        if delete:
            if len(columns) > matrix_obj.dims[1]:
                raise IOError('More columns to be deleted %d than exist %d in input file %s' % (len(columns), matrix_obj.dims[1], source))
            
            new_data = numpy.zeros((matrix_obj.dims[0], matrix_obj.dims[1]-len(columns)), dtype=numpy.double)
            new_labels = []
            new_units  = []

            new_col_idx = 0
            for old_col_idx in range(matrix_obj.dims[1]):
                if old_col_idx not in columns:
                    new_labels.append(matrix_obj.labels[old_col_idx])
                    new_units.append(matrix_obj.units[old_col_idx])

                    new_data[:,new_col_idx] = matrix_obj.data[:,old_col_idx]

                    new_col_idx += 1

            matrix_obj.data = new_data
            matrix_obj.labels = new_labels
            matrix_obj.units  = new_units

        if modify != None and len(modify) > 0:
            
            modifyDict = copy_module.copy(valuesDict)

            Get_Constant_Values(modifySect.Get_Section('->CONSTANTS'), modifyDict)

            for row_idx in rows:
                for col_idx in columns:
                    modifyDict['original'] = str(matrix_obj.data[row_idx][col_idx])

                    modify_str = Apply_Template(modify, modifyDict, mapDict=mapDict)

                    try:
                        matrix_obj.data[row_idx][col_idx] = eval(modify_str)
                    except:
                        raise RuntimeError('Error evaluating modify string: "%s"' % modify_str)

    matrix_obj.write(destination, auto_size_cols=False)
    
