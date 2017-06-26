import os
import sys
import logging
from types import ListType

from OCO_Matrix import OCO_Matrix
from OCO_TextUtils import evaluate_bool_str
from Generator_Utils import *

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict, buffer_objs):
    logger = logging.getLogger(os.path.basename(__file__))

    logger.debug('Creating aerosol blocks for file: %s' % source)

    fileObj = L2_Input.Input_File(source)

    param_def_sec = fileObj.Get_Section('PARAMETER_DEFINITION')
    if len(param_def_sec) == 0:
        logger.error('%s has sections: %s' % (source, fileObj.Get_All_Section_Names()))
        raise IOError('Could not find PARAMETER_DEFINITION section in source: %s' % source)

    orig_aero_secs     = param_def_sec[0].Get_Section('->AEROSOL')
    orig_aero_sec_names = [ curr_aer_sec.Get_Keyword_Value('name').upper() for curr_aer_sec in orig_aero_secs ]

    try:
        aero_block_prototype = orig_aero_secs[0]
    except IndexError:
        # For now just raise the error
        raise IndexError('Could not find aerosol block to use as prototype in file: %s' % source)
        
        logger.warning('No aerosol block found as prototype, trying to create a new section')
       
        aero_block_prototype = L2_Input.Section('AEROSOL')
        aero_block_prototype.Set_Keyword_Value('name', None)
        aero_block_prototype.Set_Keyword_Value('a_priori', None)
        aero_block_prototype.Set_Keyword_Value('covariance', None)
        aero_block_prototype.Set_Keyword_Value('mie_file', None)
        aero_block_prototype.Set_Keyword_Value('moment_file', None)
        aero_block_prototype.Set_Keyword_Value('retrieval_indicies', None)
    
    for curr_sect_index, curr_section_obj in enumerate(moduleSections):
        profile_names = Apply_Template(curr_section_obj.Get_Keyword_Value('profile_names'), valuesDict, mapDict=mapDict)
        type_names    = Apply_Template(curr_section_obj.Get_Keyword_Value('type_names'), valuesDict, mapDict=mapDict)

        aerosol_file  = Apply_Template(curr_section_obj.Get_Keyword_Value('from_file'), valuesDict, mapDict=mapDict)
        set_retrieval_vector = evaluate_bool_str( curr_section_obj.Get_Keyword_Value('set_retrieval_vector') )

        if curr_sect_index == 0:
            remove_existing_blocks = True
        else:
            remove_existing_blocks = False

        if aerosol_file != None and (profile_names == None or len(profile_names) == 0):
            if buffer_objs.has_key(aerosol_file):
                aerosol_obj = buffer_objs[aerosol_file]
                aerosol_obj.seek(0)
            else:
                aerosol_obj = aerosol_file

            profile_names = []
            type_names = []

            logger.debug('Using aerosol file for profile names: %s' % aerosol_obj)

            mat_obj = OCO_Matrix(aerosol_obj)
            for lbl_name in mat_obj.labels_lower:
                if lbl_name != 'pressure':
                    profile_names.append(lbl_name.upper())
                    type_names.append(lbl_name.lower())
        else:
            if profile_names == None or len(profile_names) == 0:
                raise AttributeError('profile_names needs to be defined')
            elif not type(profile_names) is ListType:
                profile_names = profile_names.split()

            if type_names == None or len(type_names) == 0:
                raise AttributeError('type_names needs to be defined')
            elif not type(type_names) is ListType:
                type_names = type_names.split()

        logger.debug('Using profile names: %s' % (', '.join(profile_names)))
        logger.debug('Using type names: %s' % (', '.join(type_names)))

        if remove_existing_blocks:
            param_def_sec[0].Set_Keyword_Value('aerosol_types', profile_names)
        else:
            existing_types = param_def_sec[0].Get_Keyword_Value('aerosol_types')

            if hasattr(existing_types, '__iter__'):
                existing_types += profile_names
            else:
                existing_types += ' ' + ' '.join(profile_names)

            param_def_sec[0].Set_Keyword_Value('aerosol_types', existing_types)

        if set_retrieval_vector:
            retrieval_vector = param_def_sec[0].Get_Keyword_Value('retrieval_vector')

            if retrieval_vector == None:
                raise IOError('Could not find retrieval_vector keyword for PARAMETER_DEFINITION in file: %s' % source)
            
            if type(retrieval_vector) is not ListType:
                retrieval_vector = retrieval_vector.split()

            retrieval_vector = [ rv_str.upper() for rv_str in retrieval_vector ]

            # Remove existing aerosol names from retrieval vector
            for curr_aer_name in orig_aero_sec_names:
                if curr_aer_name in retrieval_vector:
                    retrieval_vector.remove(curr_aer_name)

            for curr_prof_name in profile_names:
                retrieval_vector.append(curr_prof_name)

            param_def_sec[0].Set_Keyword_Value('retrieval_vector', retrieval_vector)

        # Delete from param list aerosol types
        if remove_existing_blocks:
            for as_del in orig_aero_secs:
                param_def_sec[0].children.remove(as_del)

        for (curr_prof_name, curr_type_name) in zip(profile_names, type_names):
            new_aero = copy.deepcopy(aero_block_prototype)

            new_aero.Set_Keyword_Value('name', curr_prof_name)

            mie = new_aero.Get_Keyword_Value('mie_file')
            mie_dn = os.path.dirname(mie)
            new_aero.Set_Keyword_Value('mie_file', mie_dn + '/' + curr_type_name + '.mie')

            mom = new_aero.Get_Keyword_Value('moment_file')
            mom_dn = os.path.dirname(mom)
            new_aero.Set_Keyword_Value('moment_file', mom_dn + '/' + curr_type_name + '.mom')

            param_def_sec[0].children.append(new_aero)

    logger.debug('Writing aerosol changes to file: %s' % destination)
    fileObj.Write(destination)
