import os
import sys
import logging

from Generator_Utils import NamedStringIO

from Generator_Utils import *

from OCO_TextUtils import evaluate_bool_str

from OCO_Matrix import OCO_Matrix

# Use capability of script in L2_Support/utils
from noisify_spectra import noisify_spectra_obj

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict, buffer_objs):

    matrix_obj = OCO_Matrix(source)
    for noisifySect in moduleSections:       
        noise_cut_off = Apply_Template(noisifySect.Get_Keyword_Value('noise_cut_off'), valuesDict, mapDict=mapDict)
        pixel_rows    = Apply_Template(noisifySect.Get_Keyword_Value('pixel_rows'), valuesDict, mapDict=mapDict)

        noisify_spectra_obj(matrix_obj, row_range_spec=pixel_rows, noise_cut_off=noise_cut_off)
        
    matrix_obj.write(destination, auto_size_cols=False)
