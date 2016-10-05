import os
import sys
import logging

from types import ListType 

from Generator_Utils import *

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    logger = logging.getLogger(os.path.basename(__file__))
    
    if len(moduleSections) > 1:
        raise RuntimeError('Only one CONTENTS per FILE')

    if type(destination) is str:
        if os.path.exists(destination):
            raise IOError('Source file: %s already exists, will not overwrite. Use PICK module instead' % destination)

        out_obj = open(destination, 'w')
    elif hasattr(destination, 'write'):
        out_obj = destination
    else:
        raise Exception('Unhandled object type: %s' % out_obj)

    content_lines = moduleSections[0].Get_Matrix_Data()

    logger.debug('Writing contents to %s' % destination)

    for curr_line in content_lines:
        if type(curr_line) is ListType:
            curr_line = ' '.join(curr_line)

        curr_line = Apply_Template(curr_line, valuesDict, mapDict=mapDict)
        print >>out_obj, curr_line

    if type(destination) is str:
        out_obj.close()
