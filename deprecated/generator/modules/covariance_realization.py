#!/usr/bin/env python

import os
import math
import sys
import random
import logging

from types import ListType

from Generator_Utils import *
from OCO_TextUtils import index_range_list
from OCO_Matrix import OCO_Matrix

from numpy import linalg, dot, sqrt

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    logger = logging.getLogger(os.path.basename(__file__))

    # Load existing file
    matrix_obj = OCO_Matrix(source)

    for realizeSect in moduleSections:

        # Add ability to specify cols individually or using a * to goto end
        covariance = Apply_Template(realizeSect.Get_Keyword_Value('covariance'), valuesDict, mapDict=mapDict)
        column     = Apply_Template(realizeSect.Get_Keyword_Value('column'), valuesDict, mapDict=mapDict)

        if type(column) is ListType:
            raise TypeError('Only one column can be modified per file')

        if covariance == None or len(covariance) == 0:
            raise IOError('covariance file is not specified')

        if not os.path.exists(covariance):
            raise IOError('covariance file does not exist: %s' % covariance)

        cov_obj = OCO_Matrix(covariance)

        rand_factors = [ random.normalvariate(0,1) for i in range(matrix_obj.dims[0]) ]
        (eigen_val, eigen_vec) = linalg.eigh(cov_obj.data)

        try:
            column_idx = int(column)
        except:
            if column == None:
                raise IOError('column named not defined for source file: %s' % (source))
            elif not column.lower() in matrix_obj.labels_lower:
                raise IOError('column named %s not found in source file: %s' % (column, source))

            column_idx = matrix_obj.labels_lower.index(column.lower())

        update_vals = dot(eigen_vec.transpose(), (rand_factors* sqrt(eigen_val)))
        matrix_obj.data[:, column_idx] = matrix_obj.data[:, column_idx] + update_vals

    matrix_obj.write(destination)
