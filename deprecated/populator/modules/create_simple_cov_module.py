import os
import sys
import logging

from Generator_Utils import *

from OCO_TextUtils import evaluate_bool_str

# Use capability of script in L2_Support/utils
from create_simple_cov import create_simple_cov

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    logger = logging.getLogger(os.path.basename(__file__))

    if len(moduleSections) > 1:
        raise RuntimeError('Only one section allowed per file')

    scaling = Apply_Template(moduleSections[0].Get_Keyword_Value('scaling'), valuesDict, mapDict=mapDict)
    columns = Apply_Template(moduleSections[0].Get_Keyword_Value('columns'), valuesDict, mapDict=mapDict)

    if scaling == None:
        raise ValueError('scaling keyword not defined')

    if columns == None:
        raise ValueError('cols keyword not defined')
   
    create_simple_cov(source, destination, scaling, columns)
