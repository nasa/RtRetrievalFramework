import os
import sys
import logging

from Generator_Utils import *

from OCO_TextUtils import evaluate_bool_str

# Use capability of script in L2_Support/utils
from resample_aerosol import resample_aerosol_file

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):

   if len(moduleSections) > 1:
       raise RuntimeError('Only one resample aerosol block allowed')

   resample_to = Apply_Template(moduleSections[0].Get_Keyword_Value('resample_to'), valuesDict, mapDict=mapDict)

   if resample_to == None or len(resample_to) == 0:
      raise ValueError('resample_to keyword not specified')

   if resample_to.isdigit():
      raise ValueError('resample_to must be a pressure filename not a number of levels')

   if not os.path.exists(resample_to):
      raise IOError('could not find destination pressure grid file: %s' % resample_to)

   resample_aerosol_file(source, destination, resample_to)
