import os

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix
from Generator_Utils import *

from match_fts_oco_times import match_fts_oco_times

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):

    if len(moduleSections) > 1:
        raise RuntimeError('Only one input file config set per file')
    
    if str(source) == str(destination):
        raise ValueError('source and dest filenames must be different. will not overwrite source file')

    sounding_id_list = Apply_Template(moduleSections[0].Get_Keyword_Value('sounding_id_list'), valuesDict, mapDict=mapDict)
    runlog_file      = Apply_Template(moduleSections[0].Get_Keyword_Value('runlog_file'), valuesDict, mapDict=mapDict)

    match_fts_oco_times(sounding_id_list, runlog_file, destination)
