import os
import logging

from OCO_Matrix import OCO_Matrix
from OCO_TextUtils import evaluate_bool_str
from Generator_Utils import *
import ACOS_File


from extract_gosat_noise_from_hdf import extract_sounding_noise, CNV_COLUMN_TMPL, NOISE_VALUES_HEADER_NAME

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    logger = logging.getLogger(os.path.basename(__file__))

    for curr_section_obj in moduleSections:
        hdf_filename = Apply_Template(curr_section_obj.Get_Keyword_Value('hdf_filename'), valuesDict, mapDict=mapDict)
        sounding_id = Apply_Template(curr_section_obj.Get_Keyword_Value('sounding_id'), valuesDict, mapDict=mapDict)

        # Load up parameters for modification
        param_sections = curr_section_obj.Get_Section('->PARAMS')

        # Load up an ACOS_File.L1B
        l1b_obj = ACOS_File.L1B(hdf_filename)

        if param_sections != None and len(param_sections) > 0:
            param_strs = param_sections[0].Get_Matrix_Data()
            parameters = []
            for band_strs in param_strs:
                parameters.append( numpy.array([ float(curr_val) for curr_val in band_strs ]) )
        else:
            raise Exception('PARAMS section missing')
            parameters = None

        # Extract noise from hdf file into ascii file
        try:
            dest_obj = extract_sounding_noise(l1b_obj, sounding_id, modification_params=parameters)
            logger.debug('Writing modified GOSAT noise file: %s' % destination)
            dest_obj.write(destination)
        except ValueError as val_err:
            logger.error("Failed to modify noise for sounding: %s see next message:" % sounding_id)
            logger.error(val_err)
            logger.error("Code will fail on this sounding!")
