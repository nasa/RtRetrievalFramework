import os
import sys
import logging

from Generator_Utils import *

from OCO_TextUtils import evaluate_bool_str

# Use capability of script in L2_Support/utils
from interpol_cov import interpol_cov

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):

    if len(moduleSections) > 1:
        raise RuntimeError('Only one section allowed per file')

    src_pressure_file = Apply_Template(moduleSections[0].Get_Keyword_Value('src_pressure_file'), valuesDict, mapDict=mapDict)
    dst_pressure_file = Apply_Template(moduleSections[0].Get_Keyword_Value('dst_pressure_file'), valuesDict, mapDict=mapDict)

    if src_pressure_file == None:
        raise ValueError('src_pressure_file keyword not defined')

    if dst_pressure_file == None:
        raise ValueError('dst_pressure_file keyword not defined')
   
    interpol_cov(source, destination, src_pressure_file, dst_pressure_file)
