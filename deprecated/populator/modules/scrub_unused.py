import os
import sys
import logging
from types import ListType

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix
from Generator_Utils import *

from OCO_TextUtils import evaluate_bool_str

run_basename = 'oco_l2.run'

def get_obj_referenced_filenames(section_obj, obj_dir=None):

    if obj_dir != None:
        orig_dir = os.getcwd()
        os.chdir(obj_dir)

    found_filenames=[]
    for child in section_obj.children:
        if child.type == 'assignment':
            child_val = section_obj.Get_Keyword_Value(child.leaf)

            if not type(child_val) is ListType:
                chk_filename = os.path.realpath(child_val)
                
                if os.path.exists(chk_filename):
                    found_filenames.append(chk_filename)
        elif child.type == 'section':
            for sect_filename in get_obj_referenced_filenames(child, None):
                found_filenames.append(sect_filename)

    if obj_dir != None:
        os.chdir(orig_dir)

    return found_filenames

def get_all_filenames(base_dir):

    found_files = []
    for dir_file in os.listdir(base_dir):
        real_file = os.path.realpath(os.path.join(base_dir, dir_file))
        if not os.path.islink(real_file) and os.path.isdir(real_file):
            for sub_file in get_all_filenames(real_file):
                found_files.append(sub_file)
        else:
            found_files.append(real_file)

    return found_files

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    logger = logging.getLogger(os.path.basename(__file__))
    
    if len(moduleSections) > 1:
        raise RuntimeError('Only one of this module per FILE')

    verbose  = evaluate_bool_str( moduleSections[0].Get_Keyword_Value('verbose'), False)
    max_dirs = moduleSections[0].Get_Keyword_Value('max_dirs')

    if max_dirs != None:
        max_dirs = int(max_dirs)

    run_spec_obj = L2_Input.Input_File(source)
    run_dirs = run_spec_obj.Get_Matrix_Data()

    # No directories to search for scrubbing
    if run_dirs == None:
        logger.info('No run directories available for searching for scrubbing')
        return

    sys.stdout.write('%s: Searching run directories for unreferenced files: ' % os.path.basename(__file__))

    ref_unique_filenames = {}
    dir_count = 0
    for curr_dir in run_dirs:
        if max_dirs != None and dir_count >= max_dirs:
            break
                
        run_filename = '%s/%s' % (curr_dir, run_basename)
        run_file_obj = L2_Input.Input_File(run_filename)
        control_sect = run_file_obj.Get_Section('CONTROL')
        if len(control_sect) == None:
            raise IOError('No CONTROL section in %s' % run_filename)

        try:
            input_file   = '%s/%s' % (curr_dir, control_sect[0].Get_Keyword_Value('input_file'))
        except:
            raise LookupError('Could not find find input_file keyword in CONTROL section of file: %s' % run_filename)
        inp_file_obj = L2_Input.Input_File(input_file)

        run_file_list = get_obj_referenced_filenames(run_file_obj, curr_dir)
        inp_file_list = get_obj_referenced_filenames(inp_file_obj, curr_dir)

        for ref_filename in (run_file_list + inp_file_list):
            if ref_unique_filenames.has_key(ref_filename):
                ref_unique_filenames[ref_filename] += 1
            else:
                ref_unique_filenames[ref_filename] = 0

        # Progress marks since this loop takes awhile
        sys.stdout.write('.')
        sys.stdout.flush()
        
        dir_count += 1

    # Seperate progress marks
    sys.stdout.write('\n')

    if type(destination) is str:
        dest_filenames = get_all_filenames(destination)
    else:
        dest_filenames = get_all_filenames(destination.filename)

    for static_filename in dest_filenames:
        if not ref_unique_filenames.has_key(static_filename):
            try:
                os.remove(static_filename)
                if verbose:
                    logger.debug('Deleted: %s' % static_filename)
            except:
                logger.error('Could not delete: %s' % static_filename)
