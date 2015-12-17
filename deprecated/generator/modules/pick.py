import copy
import logging

import L2_Input
from types import ListType

from Generator_Utils import *
from OCO_TextUtils import index_range_list, evaluate_bool_str

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    logger = logging.getLogger(os.path.basename(__file__))
    
    logger.debug('')
    logger.debug('Reading source file: %s' % source)
    modFileObj = L2_Input.Input_File(source)
    
    for pick in moduleSections:
        section  = Apply_Template(pick.Get_Keyword_Value('section'), valuesDict, mapDict=mapDict)
        keyword  = Apply_Template(pick.Get_Keyword_Value('keyword'), valuesDict, mapDict=mapDict)

        # Resolve template later when can add to the values dictionary
        template = pick.Get_Keyword_Value('template', addWhiteSpace=True)

        which_line    = Apply_Template(pick.Get_Keyword_Value('which_line'), valuesDict, mapDict=mapDict)
        which_section = Apply_Template(pick.Get_Keyword_Value('which_section'), valuesDict, mapDict=mapDict)
        which_keyword = Apply_Template(pick.Get_Keyword_Value('which_keyword'), valuesDict, mapDict=mapDict)

        ignore_missing = Apply_Template(pick.Get_Keyword_Value('ignore_missing'), valuesDict, mapDict=mapDict)
        index_format   = Apply_Template(pick.Get_Keyword_Value('index_format'), valuesDict, mapDict=mapDict)

        delete   = Apply_Template(pick.Get_Keyword_Value('delete'), valuesDict, mapDict=mapDict)
        indent   = Apply_Template(pick.Get_Keyword_Value('indent'), valuesDict, mapDict=mapDict)
        unique   = evaluate_bool_str(Apply_Template(pick.Get_Keyword_Value('unique'), valuesDict, mapDict=mapDict), False)

        if indent == None:
            indent = 0

        if keyword != None and type(keyword) is not ListType:
            keyword = [keyword]

        if which_line == None:
            which_line = None
        elif type(which_line) is ListType or not which_line.isdigit():
            raise ValueError('which_line must be a scalar integer')
        else:
            which_line = int(which_line)

        if (template == None or len(template) == 0) and (delete == None or len(delete)) == 0:
            raise ValueError('template must be defined for PICK')

        if section != None:
            if keyword != None:
                if delete != None and len(delete) > 0:
                    logger.debug('Deleting keyword %s->%s' % (section, keyword))
                else:
                    logger.debug('Modifying keyword %s->%s' % (section, keyword))
            elif delete != None and len(delete) > 0:
                logger.debug('Deleting section %s' % section)
        else:
            if keyword != None:
                if delete != None and len(delete) > 0:
                    logger.debug('Deleting keyword %s' % keyword)
                else:
                    logger.debug('Modifying keyword %s' % keyword)

            elif delete != None and len(delete) > 0:
                logger.debug('Deleting lines from root file section:', delete)

        # Find the section to modify
        if section == None:
            modSect = [ modFileObj.rootNode ]
        else:
            if type(section) is ListType:
                modSect = []
                for curr_sect_name in section:
                    for found_sect in modFileObj.rootNode.Get_Section(curr_sect_name):
                        modSect.append(found_sect)
            else:
                modSect = modFileObj.rootNode.Get_Section(section)

        if len(modSect) == 0:
            modSect = [L2_Input.Section(leaf=section)]

            if which_line != None:
                modFileObj.children.insert(which_line, modSect[0])
            else:
                modFileObj.children.append(modSect[0])

        # If which is defined then filter sections to modify
        if keyword != None and which_section != None:

            try:
                section_indexes = index_range_list(which_section)
            except:
                section_indexes = []
                curr_index = 0
                for testSect in modSect:
                    sectName = testSect.Get_Keyword_Value('name')
                    if sectName == which_section:
                        section_indexes.append(curr_index)
                    curr_index += 1
                if len(section_indexes) == 0 and not ignore_missing:
                    raise IOError('Could not find section named: %s with name keyword: %s in file %s' % (section, which_section, source))
                
            sectChoices = []
            for w_idx in section_indexes:
                try:
                    sectChoices.append( modSect[w_idx] )
                except:
                    raise IOError("Section index: %d not found for keyword: %s, section: %s in file %s" % (w_idx, keyword, section, source))
            modSect = sectChoices

        # Finally modify all chosen sections and chosen keywords
        for pickSect in modSect:
            pickValDict = copy.copy(valuesDict)

            if keyword != None:
                for curr_keyname in keyword:
                    keyword_val = pickSect.Get_Keyword_Value(curr_keyname)
                    if type(keyword_val) is ListType:
                        pickValDict[curr_keyname] = ' '.join(keyword_val)
                    elif keyword_val != None:
                        pickValDict[curr_keyname] = keyword_val

            newValue = Apply_Template(template, pickValDict, mapDict=mapDict)

            if delete != None and len(delete) > 0:

                try:
                    delete = index_range_list(delete)
                except:
                    delete = [0]

                if keyword != None and len(keyword) > 0:
                    for curr_keyname in keyword:
                        pickSect.Delete_Keyword(curr_keyname, which=delete)
                elif section != None and len(section) > 0:
                    modFileObj.Delete_Section(pickSect) # which=delete
                else:
                    # Remove from a location a certain # of times
                    for dlist_idx in range(len(delete)-1, -1, -1):
                        del_idx = delete[dlist_idx]
                        x = pickSect.children.pop(del_idx)

            elif keyword != None:

                if which_keyword != None:
                    which_keyword = index_range_list(which_keyword)

                for curr_keyname in keyword:
                    if index_format != None and hasattr(newValue, '__iter__'):
                        for curr_index, curr_value in enumerate(newValue):
                            index_keyname = index_format.format(keyword=curr_keyname, index=curr_index+1)
                            pickSect.Set_Keyword_Value(index_keyname, curr_value, which=which_keyword, indent=indent)
                    else:
                        pickSect.Set_Keyword_Value(curr_keyname, newValue, which=which_keyword, indent=indent)
            else:
                if not type(newValue) is ListType:
                    newValue = [ newValue ]

                for currValue in newValue:
                    if unique and currValue in modFileObj.Get_Matrix_Data():
                        continue
                    
                    newNode = L2_Input.Node('value', currValue + '\n')

                    if which_line != None:
                        pickSect.children.insert(which_line, newNode)
                    else:
                        pickSect.children.append(newNode)

    # Write newly modified file
    logger.debug('Writing destination file: %s' % destination)
    modFileObj.Write(destination)
