import os
import re
import sys
import logging
from types import ListType

import numpy
from OCO_Matrix import OCO_Matrix
from OCO_TextUtils import evaluate_bool_str
from Generator_Utils import *

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    logger = logging.getLogger(os.path.basename(__file__))
    
    if len(moduleSections) > 1:
        raise AttributeError('only one merge_profile allowed per file')
    
    file_obj = OCO_Matrix(source)

    # Make sure not adding logarithmic aerosol data
    in_log = False
    header_keys = [ head_str.lower() for head_str in file_obj.header.keys() ]
    if 'retrieval_mode' in header_keys and file_obj.header[file_obj.header.keys()[header_keys.index('retrieval_mode')]] == 'logarithmic':
        in_log = True
    
    src_data           = file_obj.data
    src_profile_names  = file_obj.labels_lower
    src_units          = file_obj.units

    src_pressure_col   = src_profile_names.index('pressure')
        
    dst_profile_names = moduleSections[0].Get_All_Keyword_Names()

    dst_labels = {}
    dst_units  = {}
    dst_data   = {}

    for dest_profile_name in dst_profile_names:
        source_profiles_srch = Apply_Template(moduleSections[0].Get_Keyword_Value(dest_profile_name), valuesDict, mapDict=mapDict)

        if type(source_profiles_srch) is not ListType:
            source_profiles_srch = [ source_profiles_srch ]

        for source_profile_re in source_profiles_srch:
            for src_profile in src_profile_names:
                if re.search(source_profile_re, src_profile):
                    src_col_idx = src_profile_names.index(src_profile)

                    # Don't ever use the pressure column even if it would match
                    if src_col_idx == src_pressure_col:
                        continue

                    if not dst_data.has_key(dest_profile_name):
                        dst_data[dest_profile_name]   = numpy.zeros(src_data.shape[0], dtype=float)
                        dst_units[dest_profile_name]  = None
                        dst_labels[dest_profile_name] = dest_profile_name

                    if in_log:
                        dst_data[dest_profile_name] += numpy.exp(src_data[:, src_col_idx])
                    else:
                        dst_data[dest_profile_name] += src_data[:, src_col_idx]

                    # Complain if mixing units when adding
                    if src_col_idx < len(src_units):
                        if dst_units[dest_profile_name] == None:
                            dst_units[dest_profile_name] = src_units[src_col_idx]
                        elif dst_units[dest_profile_name].lower() != src_units[src_col_idx].lower():
                            raise Exception('Mixing unit types when merging src profile %s into dst profile %s' % (src_profile, dest_profile_name))

    if 'pressure' in file_obj.labels_lower:
        file_obj.labels = ['pressure']
        file_obj.units  = ['Pa']
        new_data = numpy.zeros((src_data.shape[0], len(dst_data.keys())+1), dtype=float)
        new_data[:, 0] = file_obj.data[:, src_pressure_col]
        beg_col = 1
    else:
        file_obj.labels = []
        file_obj.units  = []
        new_data = numpy.zeros((src_data.shape[0], len(dst_data)), dtype=float)
        beg_col = 0

    for dest_profile_name, dest_col_index in zip(dst_data.keys(), range(beg_col, beg_col+len(dst_data))):
        if in_log:
            new_data[:, dest_col_index] = numpy.log(dst_data[dest_profile_name][:])
        else:
            new_data[:, dest_col_index] = dst_data[dest_profile_name][:]
        file_obj.labels.append( dst_labels[dest_profile_name] )
        file_obj.units.append( dst_units[dest_profile_name] )

    logger.debug('Merged profile names: %s' % file_obj.labels[1:])

    file_obj.data = new_data
    file_obj.write(destination)
