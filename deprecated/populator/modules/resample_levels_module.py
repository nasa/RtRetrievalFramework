import os
import sys
import logging

from Generator_Utils import *

from OCO_TextUtils import evaluate_bool_str

# Use capability of script in L2_Support/utils
from resample_levels import resample_levels

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):

   if len(moduleSections) > 1:
       raise RuntimeError('Only one resample levels block allowed')

   resample_to = Apply_Template(moduleSections[0].Get_Keyword_Value('resample_to'), valuesDict, mapDict=mapDict)
   extrapolate = evaluate_bool_str(moduleSections[0].Get_Keyword_Value('extrapolate'))

   if resample_to == None or len(resample_to) == 0:
       raise ValueError('resample_to keyword not specified')

   resample_levels(source, destination, resample_to, val_extrapolate=extrapolate)
