import re
import logging
from types import ListType

import L2_Input
from Generator_Utils import *

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    logger = logging.getLogger(os.path.basename(__file__))

    fileObj = L2_Input.Input_File(source)
    
    for substSect in moduleSections:
        from_spec = Apply_Template(substSect.Get_Keyword_Value('from'), valuesDict, mapDict=mapDict)
        to_spec   = Apply_Template(substSect.Get_Keyword_Value('to'), valuesDict, mapDict=mapDict)
        not_str   = Apply_Template(substSect.Get_Keyword_Value('not'), valuesDict, mapDict=mapDict)
        keyword   = Apply_Template(substSect.Get_Keyword_Value('keyword'), valuesDict, mapDict=mapDict)
        section   = Apply_Template(substSect.Get_Keyword_Value('section'), valuesDict, mapDict=mapDict)

        if to_spec == None:
            to_spec = ''

        if from_spec == None:
            raise ValueError('from keyword must be defined')

        if section != None and type(section) is not ListType and section.find('->') >= 0:
            search_sect = fileObj.Get_Section(section)

            if len(search_sect) == 0:
                raise IOError('Could not find section %s in file: %s' % (section, source))
        else:
            search_sect = [ fileObj.rootNode ]

        if type(from_spec) is not ListType:
            from_list = [ from_spec ]
        else:
            from_list = from_spec

        for from_str in from_list:
            for curr_sect in search_sect:
                Substitute_Node(curr_sect, from_str, to_spec, not_str, keyword, section)
      
    fileObj.Write(destination)


def Substitute_Node(node, from_str, to_str, not_str=None, keyword=None, section=None, match_section=False):

    if node.type == 'section':

        match_section = False
        if section == None:
            match_section = True
        else:
            if type(section) is not ListType:
                sect_list = [ section ]
            else:
                sect_list = section

            for curr_sect in sect_list:
                if re.search(curr_sect, node.leaf[0]):
                    match_section = True
            
        for child in node.children:
            Substitute_Node(child, from_str, to_str, not_str, keyword, section, match_section)

    elif node.type == 'assignment' and match_section:

        match_keyword = False
        if keyword == None:
            match_keyword = True
        else:
            if type(keyword) is not ListType:
                keyword_list = [ keyword ]
            else:
                keyword_list = keyword

            for curr_key in keyword_list:
                if re.search(curr_key, node.leaf):
                    match_keyword = True
                    break

        if match_keyword: 
            for child in node.children:
                if child.leaf != None and not_str == None or not re.search(not_str, child.leaf):

                    match_obj = re.search(from_str, child.leaf)
                    if match_obj:
                        beg = match_obj.start()
                        end = match_obj.end()

                        match_from_str = child.leaf[beg:end]
                        child.leaf = child.leaf.replace(match_from_str, to_str)
