import os
import logging

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix
from Generator_Utils import *

pout_col_idx   = 27
convert_factor = 1e2

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):

    if len(moduleSections) > 1:
        raise RuntimeError('Only one input file config set per file')
    
    if str(source) == str(destination):
        raise ValueError('source and dest filenames must be different. will not overwrite source file')

    spectrum_file = Apply_Template(moduleSections[0].Get_Keyword_Value('spectrum_file'), valuesDict, mapDict=mapDict)


    if type(source) is str and not os.path.exists(source):
        raise IOError('Runlog file %s does not exist' % source)

    base_spec_name = os.path.basename(spectrum_file)

    # Use grep because its faster than doing it outself
    matched_line = None
    if type(source) == str:
        grep_cmd = "grep -E " + base_spec_name + " " + source
        matched_line = os.popen(grep_cmd).readline()

    elif hasattr(source, 'read'):
        for curr_line in source.readlines():
            if re.search(base_spec_name, curr_line):
                matched_line = curr_line
                break
    else:
        raise Exception('Unsupported object: %s' % source)

    if matched_line == None or len(matched_line) == 0:
        raise IOError('Could not find spectrum name: %s in run log file: %s' % (base_spec_name, source))

    try:
        matched_columns = matched_line.split()
        psurf_val = float(matched_columns[pout_col_idx]) * convert_factor
    except:
        raise ValueError('Failed to parse psurf value from: "%s" from runlog line: %s' % (matched_columns[pout_col_idx], matched_line))
   
    out_obj = OCO_Matrix()
    
    out_obj.data = numpy.zeros((1,1), dtype=float)
    out_obj.data[0,0] = psurf_val

    out_obj.file_id = 'psurf value extracted for spectrum named: %s from runlog file: %s' % (base_spec_name, source)
    out_obj.labels = ['PSURF']

    out_obj.write(destination)
