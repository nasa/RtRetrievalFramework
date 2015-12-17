import os
import sys
import logging

from Generator_Utils import *

from OCO_Matrix import OCO_Matrix
from OCO_TextUtils import evaluate_bool_str

from H5Dump_Util import *
import Orbit_Sim

# Use capability of script in L2_Support/utils
from extract_orbit_sim_data import extract_orbit_sim_data
  
def Process_File(sourceFilename, destFilename, fileKeywords, moduleSections, valuesDict, mapDict):
   logger = logging.getLogger(os.path.basename(__file__))

   if len(moduleSections) > 1:
       raise RuntimeError('Only one extraction block allowed')

   prof_file    = Apply_Template(moduleSections[0].Get_Keyword_Value('prof_file'), valuesDict, mapDict=mapDict)
   log_file     = Apply_Template(moduleSections[0].Get_Keyword_Value('log_file'), valuesDict, mapDict=mapDict)
   l1b_file     = Apply_Template(moduleSections[0].Get_Keyword_Value('l1b_file'), valuesDict, mapDict=mapDict)
   resample_to  = Apply_Template(moduleSections[0].Get_Keyword_Value('resample_to'), valuesDict, mapDict=mapDict)
   verbose      = evaluate_bool_str(moduleSections[0].Get_Keyword_Value('verbose'))
   id_list_file = Apply_Template(moduleSections[0].Get_Keyword_Value('id_list_file'), valuesDict, mapDict=mapDict)
   id_section   = Apply_Template(moduleSections[0].Get_Keyword_Value('id_section'), valuesDict, mapDict=mapDict)

   if not os.path.isdir(destFilename):
      raise IOError('destFilename %s must be a directory not a file' % destFilename)

   if id_list_file != None and len(id_list_file) > 0:
      sounding_ids = [ long(id_val) for id_val in Read_Id_List_File(id_list_file, id_section) ]
   else:
      sounding_ids = None

   if l1b_file == None:
      raise ValueError('l1b_file is must be defined')

   if log_file != None:
      logger.debug('Extracting orbit simulator data from %s into %s' % (prof_file, destFilename))
   else:
      logger.debug('Extracting orbit simulator data from %s and %s into %s' % (prof_file, log_file, destFilename))
   
   extract_orbit_sim_data(prof_file, l1b_file, destFilename, log_file=log_file, resample_to=resample_to, sounding_id_list=sounding_ids, verbose=verbose)
