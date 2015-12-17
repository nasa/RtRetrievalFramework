import os
import re
import sys
import copy
import glob
import logging
from StringIO import StringIO

# For use in eval substititions
import h5py
import numpy
import random
import math

from types import ListType, StringType

from OCO_TextUtils import index_range_list, evaluate_bool_str
import L2_Input

DEFAULT_MAP_NAME = 'default'

class TemplateError(AttributeError):
    "Exception for unresolved template substitions"

class HandledProcessingError(Exception):
    "Exception for processing errors that are handled internally and just need a message passed up"

class NamedStringIO(StringIO):
    """Gives a name to a StringIO. Does not subclass so it can use cStringIO for speed"""

    def __init__(self, filename=None):
        if filename != None:
            self.filename = filename
        else:
            self.filename = None

        if self.filename != None:
            StringIO.__init__(self, open(self.filename).read())
        else:
            StringIO.__init__(self)

    def __eq__(self, other):
        if other == None:
            return False
        else:
            return self.filename == other.filename

    def __ne__(self, other):
        return not self.__eq__(other)

    def close(self):
        raise Exception('close method should not be called by a module')

    def close_buffer(self):
        #logger = logging.getLogger('NamedStringIO::close')
        #logger.debug('Closing: %s(%s)' % (self.filename, hex(id(self))))
        StringIO.close(self)

    def write_file(self, filename=None):
        if filename != None:
            self.filename = filename
            
        self.seek(0)
        with open(self.filename, 'w') as file_obj:
            file_obj.write( self.getvalue() )

    def __str__(self):
        # Try and return a string that would fail if used as
        # a real filename, hence the the // prefix
        return '//%s(%s): "%s"//' % ("StringIO", hex(id(self)), self.filename)
        

def Apply_Template(template_spec, valuesDict, mapDict={}):

    if template_spec == None or (type(template_spec) is StringType and len(template_spec) == 0):
        return None

    if type(template_spec) is ListType:
        template_eval = []
        for sub_spec in template_spec:
            template_eval.append( Apply_Template(sub_spec, valuesDict, mapDict) )
        return template_eval
    else:
        if not type(template_spec) is StringType:
            template_spec = str(template_spec)
            
        orig_template_str = copy.copy(template_spec)
        
        # Evaluate simple substitutions
        for key_str in re.findall('<[^<>]+?>', template_spec):
            else_loc = key_str.find('|')

            key_beg = 1
            if else_loc > 1:
                else_value = key_str[else_loc+1:-1]
                key_end = else_loc
            else:
                else_value = None
                key_end = -1
                
            list_name = key_str[key_beg:key_end]

            replacement_value = None
            if not valuesDict.has_key(list_name):
                if else_value != None:
                    replacement_value = else_value
                else:
                    raise TemplateError('Could not find template replacement value for %s in string "%s"' % (key_str, orig_template_str))
            else:
                replacement_value = valuesDict[list_name]

            # Somehow a None can work its way in and should just
            # be an empty string
            if replacement_value == None:
                replacement_value = ''

            if hasattr(replacement_value, '__iter__'):
                new_spec = []
                for curr_repl in replacement_value:
                    new_spec.append(template_spec.replace(key_str, curr_repl))
                    
                if len(new_spec) == 1:
                    new_spec = new_spec[0]

                template_spec = new_spec
            else:
                template_spec = template_spec.replace(key_str, replacement_value)

    # Evaluate map substitions
    if type(template_spec) is ListType:
        template_eval = []
        for sub_spec in template_spec:
            template_eval.append( Apply_Template(sub_spec, valuesDict, mapDict) )
        return template_eval
    else:
        for key_str in re.findall('\[[a-zA-Z][\w:+]*?\]', template_spec):
            if mapDict == None or len(mapDict) == 0:
                raise TemplateError('Map substition strings found in string "%s" but no map dictionary supplied to function' % template_spec)
                
            key_contents = key_str[1:-1]
            key_parts = key_contents.split(':')
            if len(key_parts) > 1:
                (map_set, list_spec) = key_parts
            else:
                map_set = DEFAULT_MAP_NAME
                list_spec = key_contents

            if not mapDict.has_key(map_set):
                raise TemplateError('map set %s not found among map sets: [%s]' % (map_set, ', '.join(mapDict.keys())))

            num_subs_made = 0
            map_srch_key = ""
            for list_part in list_spec.split('+'):
                if valuesDict.has_key(list_part):
                    map_srch_key  += valuesDict[list_part]
                    num_subs_made += 1
                else:
                    map_srch_key += list_part

            if not mapDict[map_set].has_key(map_srch_key):
                raise TemplateError('Could not map value %s from map set %s in template string "%s"' % (map_srch_key, map_set, orig_template_str))

            map_repl_value = mapDict[map_set][map_srch_key]
            if type(map_repl_value) is ListType:
                template_spec_list = []
                for curr_repl_value in map_repl_value:
                    template_spec_list.append( template_spec.replace(key_str, curr_repl_value) )

                template_spec = template_spec_list
            else:
                template_spec = template_spec.replace(key_str, map_repl_value)

            # Evaluate any substitutions within maps
            template_spec = Apply_Template(template_spec, valuesDict, mapDict)

    # Only evaulate eval statements if previous subsititutions did not result in a list
    if type(template_spec) is StringType:
        eval_re = re.compile('eval\(.*\)')
        for eval_str in re.findall(eval_re, template_spec):
            eval_val = re.sub('^eval\(', '', eval_str)
            eval_val = re.sub('\)$', '', eval_val)

            eval_result = str( eval(eval_val) )

            template_spec = template_spec.replace(eval_str, eval_result)

    return template_spec

def Expand_Filename(filename, basePath=None):

    if filename == None:
        return [ '' ]
    
    # Expand any usage of ~/ or ~user in pathnames
    filename = os.path.expanduser(filename)
    if basePath != None:
        basePath = os.path.expanduser(basePath)
   
    if basePath != None and not filename[0] == '/' and not filename.find('..') == 0:
        fullname = basePath + '/' + filename
    else:
        fullname = filename

    matches = glob.glob(fullname)

    if len(matches) > 0:
        return matches
    else:
        return [ fullname ]

def Case_Relative_Path(filename, base_dir=None, case_dir=None):

    if case_dir == None:
        return filename

    if not os.path.exists(case_dir) and base_dir != None:
        case_under_base = os.path.join(base_dir, case_dir)
        if not os.path.exists(case_under_base):
            raise IOError('Could not find case directory as: "%s" or "%s"' % (case_dir, case_under_base))
        case_dir = case_under_base
    elif not os.path.exists(case_dir):
        raise IOError('Could not find case directory: %s' % case_dir)

    file_rel_to_case = os.path.relpath(filename, os.path.realpath(case_dir))

    if base_dir != None and filename.find(os.path.realpath(base_dir)) != 0:
        return filename
    else:
        return file_rel_to_case


def Should_Process(check_section, valuesDict, mapDict, default=True):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)
    
    do_process = default
    try:
        only_if_tmpl = check_section.Get_Keyword_Value('only_if')
        only_if_val = Apply_Template(only_if_tmpl, valuesDict, mapDict=mapDict)
    except TemplateError:
        only_if_val = False

    try:
        not_if_tmpl = check_section.Get_Keyword_Value('not_if')
        not_if_val  = Apply_Template(not_if_tmpl,  valuesDict, mapDict=mapDict)
    except TemplateError:
        not_if_val  = False

    # not_if takes precedence over only_if
    if not_if_val != None:
        if not getattr(not_if_val, '__iter__', False):
            not_if_val = [not_if_val]

        for curr_val in not_if_val:
            do_process = not evaluate_bool_str(str(curr_val), default=True)
            if not do_process:
                break

        logger.debug('Should process %s = %s : not_if string: "%s" evaluates: %s' % (check_section.leaf[0], do_process, not_if_tmpl, not_if_val))

    # If only_if_val is defined then make sure it evaluates to true
    if do_process and only_if_val != None:
        if not getattr(only_if_val, '__iter__', False):
            only_if_val = [only_if_val]

        for curr_val in only_if_val:
            do_process = evaluate_bool_str(str(curr_val), default=False)
            if not do_process:
                break

        logger.debug('Should process %s = %s : only_if string: "%s" evaluates: %s' % (check_section.leaf[0], do_process, only_if_tmpl, only_if_val))

    return do_process


def Get_Map_Values(mapSectionsList, subDict=None, existing={}):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)
    
    mapData = copy.copy(existing)

    if type(mapSectionsList) is not ListType:
        mapSectionsList = [ mapSectionsList ]

    for mapSection in mapSectionsList:                   
        mapName = mapSection.Get_Keyword_Value('name')
        if mapName == None or mapName == '':
            mapName = DEFAULT_MAP_NAME

        if not mapName in mapData.keys():
            if type(mapName) is ListType:
                raise ValueError('A named MAP section must contain map vaues in VALUES subsection')
            mapData[mapName] = {}

        required = evaluate_bool_str( mapSection.Get_Keyword_Value('required'), True )
        if subDict != None:
            mapFilename = Apply_Template( mapSection.Get_Keyword_Value('from_file'), subDict )
            sectionName = Apply_Template( mapSection.Get_Keyword_Value('section'), subDict )
        else:
            mapFilename = mapSection.Get_Keyword_Value('from_file')
            sectionName = mapSection.Get_Keyword_Value('section')

        if mapFilename != None and len(mapFilename) > 0:
             mapFilename = Expand_Filename(mapFilename)[0]

             if not os.path.exists(mapFilename):
                 if required:
                     raise IOError("MAP source file '%s' does not exist" % (mapFilename))
                 else:
                     continue

             if sectionName == None:
                 logger.debug('Loading MAP %s contents from file: %s' % (mapName, mapFilename))
                 mapFileData = L2_Input.Input_File(mapFilename)
                 mapSection = mapFileData.rootNode
                 mapValues = mapSection.Get_Matrix_Data()
             else:
                logger.debug('Loading MAP %s section as %s contents from file: %s' % (sectionName, mapName, mapFilename))
                fileObj = L2_Input.Input_File(mapFilename)

                foundSects = fileObj.Get_Section(sectionName)

                if foundSects == None or len(foundSects) == 0:
                    raise IOError('Could not find section %s in file: %s' % (sectionName, mapFilename))

                mapValues = []
                for currFileSect in foundSects:
                    for sectKeyName in currFileSect.Get_All_Keyword_Names():
                        sectKeyVal = currFileSect.Get_Keyword_Value(sectKeyName)
                        mapValues.append( [str(sectKeyName), str(sectKeyVal)] )

        else:
            if mapName == DEFAULT_MAP_NAME:
                valuesSect = [mapSection]
            else:
                valuesSect = mapSection.Get_Section('MAP->VALUES')

            if len(valuesSect) == 0:
                raise ValueError('Could not find map values in either main MAP section or VALUES sub section')
            if len(valuesSect) > 1:
                raise ValueError('Too many VALUES sub sections for MAP section')

            mapValues = valuesSect[0].Get_Matrix_Data()

        if mapValues != None:
            for mapRow in mapValues:
                mapKey = mapRow[0]
                mapValue = mapRow[1:]

                if len(mapValue) == 1:
                    mapValue = mapValue[0]

                mapData[ mapName ][ mapKey ] = mapValue

    return mapData

def Get_List_File_Values(listLocation, listName, sectionName=None, directoryLevels=None):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)

    # Expand any globs
    if type(listLocation) is str:
        listLocation = Expand_Filename(listLocation)

    # Make sure we have something to iterate over
    elif not hasattr(listLocation, '__iter__'):
        listLocation = [listLocation]

    if directoryLevels == None:
        directoryLevels = 1
    else:
        if not type(directoryLevels) is int and not directoryLevels.isdigit():
            raise ValueError('directoryLevels argument must be an integer or convertable to one')
        else:
            directoryLevels = int(directoryLevels)
        
    fileValues = []
    for listFile in listLocation:
       
        if type(listFile) is str and os.path.isdir(listFile):
            logger.debug('Loading LIST %s contents from directory: %s' % (listName, listFile))

            # Use directoryLevels - 1 parts of the end of the path
            old_path = listFile
            new_path = ''
            for level_idx in range(directoryLevels-1):
                old_path, dir_name = os.path.split(old_path)
                new_path = os.path.join(new_path, dir_name)

            # Use path as a list based on filenames present there
            for dir_filename in os.listdir(listFile):
                fileValues.append( os.path.join(new_path, dir_filename) ) 

        elif sectionName == None:
            logger.debug('Loading LIST %s contents from file: %s' % (listName, listFile))

            if type(listFile) is str:
                listFileObj = open(listFile, 'r')
            elif hasattr(listFile, 'read'):
                listFileObj = listFile
            else:
                raise Exception('Unknown read object: %s' % listFile)

            for listLine in listFileObj.readlines():
                if len(listLine.strip()) > 0 and listLine.strip()[0] != '#':
                    fileValues.append(listLine.strip())

            if type(listFile) is str:
                listFileObj.close()
        else:
            logger.debug('Loading LIST %s section as %s contents from file: %s' % (sectionName, listName, listFile))
            fileObj = L2_Input.Input_File(listFile)

            sectNameParts = sectionName.split('->')

            foundSects = fileObj.Get_Section('->'.join(sectNameParts[0:-1]) + '->LIST')

            fileListSect = None
            for currFileSect in foundSects:
                currListName = currFileSect.Get_Keyword_Value('name')
                if currListName != None and currListName == sectNameParts[-1]:
                    fileListSect = currFileSect.Get_Section('LIST->VALUES')
                    break

            if fileListSect == None or len(fileListSect) == 0:
                raise IOError('Could not find section %s in file: %s' % (sectionName, listFile))

            fileValues = fileListSect[0].Get_Matrix_Data()

    return fileValues

def Get_List_Section_Values(listSections, listSubDict, setMap):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)
    
    listValuesDict = {}

    if type(listSections) is not ListType:
        listSections = [ listSections ]
    
    for listNode in listSections:
        nodeName = listNode.Get_Keyword_Value('name')
        listFilename = Apply_Template( listNode.Get_Keyword_Value('from_file'), listSubDict, mapDict=setMap )
        listConstant = Apply_Template( listNode.Get_Keyword_Value('from_constant'), listSubDict, mapDict=setMap )
        sectionName  = Apply_Template( listNode.Get_Keyword_Value('section'), listSubDict, mapDict=setMap )
        directoryLevels = Apply_Template( listNode.Get_Keyword_Value('directory_levels'), listSubDict, mapDict=setMap )

        if nodeName == None or len(nodeName) == 0:
            raise ValueError('name keyword must be defined for LIST')

        valueChildren = listNode.Get_Section('LIST->VALUES')
        rangeChildren = listNode.Get_Section('LIST->RANGE')

        if listFilename != None:
            listValuesDict[nodeName] = Get_List_File_Values(listFilename, nodeName, sectionName, directoryLevels)
            
        elif listConstant != None and len(listConstant) > 0:

            listValuesDict[nodeName] = [const_val.strip() for const_val in listConstant.split()]

        elif len(rangeChildren) > 0:
            rangeSect = rangeChildren[0]
            minVal = rangeSect.Get_Keyword_Value('min')
            maxVal = rangeSect.Get_Keyword_Value('max')
            format = rangeSect.Get_Keyword_Value('format')

            if minVal == None or len(minVal) <= 0:
                minVal = 0
            else:
                minVal = int(minVal)

            if maxVal == None or len(maxVal) <= 0:
                maxVal = 0
            else:
                maxVal = int(maxVal)

            if format == None or len(format) <= 0:
                format = '%d'

            rangeList = []
            for rangeVal in xrange(minVal, maxVal):
                rangeList.append( format % rangeVal )

            listValuesDict[nodeName] = rangeList

        elif len(valueChildren) > 0:
            valueSect = valueChildren[0]
            matrixData = valueSect.Get_Matrix_Data()

            if matrixData == None:
                raise ValueError("VALUES section of LIST named '%s' can not be empty" % nodeName)

            listValuesDict[nodeName] = matrixData
        else:
            raise ValueError('LIST section named %s must contain a RANGE or VALUES subsection' % nodeName)

    # apply any substitutions to list values
    for listName in listValuesDict.keys():
        listValuesDict[listName] = Apply_Template( listValuesDict[listName], listSubDict, mapDict=setMap )

    return listValuesDict

def Get_Constant_Values(constSections, existingDict={}, mapDict={}, templateDict=None):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)

    # Use existing dict as template when template dict is not defined
    # Put for backwards compatibility, where existingDict was used
    # in keyword applications
    if templateDict == None:
        templateDict = existingDict
    
    loaded_constants = []
    for const_sect in constSections:
        sect_name = const_sect.leaf[0]
        
        for const_child in const_sect.children:
        
            if const_child.type == 'assignment':
                const_name = const_child.leaf
                const_val = Apply_Template( const_sect.Get_Keyword_Value(const_name), templateDict, mapDict=mapDict )

                if type(const_val) is ListType:
                    raise ValueError('constant %s defined more than once or is a list with value: %s' % (const_name, const_val))

                loaded_constants.append(const_name)
                existingDict[ const_name ] = const_val

            elif const_child.type == 'section' and const_child.leaf[0] == 'EXTRACT':
                extract_filename = Apply_Template( const_child.Get_Keyword_Value('filename'), templateDict, mapDict=mapDict )
                allow_missing    = evaluate_bool_str( Apply_Template( const_child.Get_Keyword_Value('allow_missing'), templateDict, mapDict=mapDict ), False)
                
                keyword_sect_list = const_child.Get_Section('EXTRACT->KEYWORDS')

                if extract_filename == None or len(extract_filename) == 0:
                    raise ValueError('filename must be specified for %s->EXTRACT section' % sect_name)

                if keyword_sect_list == None or len(keyword_sect_list) == 0:
                    raise ValueError('KEYWORD section must be specified for %s->EXTRACT section' % sect_name)

                if not os.path.exists(extract_filename):
                    raise ValueError('filename specified for %s->EXTRACT section does not exist: %s' % (sect_name, extract_filename))

                logger.debug('Reading constant keyword values from file: %s' % extract_filename)
                keyFileObj = L2_Input.Input_File(extract_filename)

                for keyword_sect in keyword_sect_list:
                    wanted_consts = keyword_sect.Get_All_Keyword_Names()
                    for new_const_name in wanted_consts:
                        search_path = Apply_Template( keyword_sect.Get_Keyword_Value(new_const_name), templateDict, mapDict=mapDict )
                        logger.debug('Loading %s from keyword file as %s' % (search_path, new_const_name))
                        search_sect_name = '->'.join(search_path.split('->')[0:-1])
                        search_key_name  = search_path.split('->')[-1]

                        search_sect_obj = keyFileObj.Get_Section(search_sect_name)

                        if search_sect_obj == None or len(search_sect_obj) == 0:
                            raise IOError('Could not find section: %s in file: %s' % (search_sect_name, extract_filename))

                        new_const_value = [ sect.Get_Keyword_Value(search_key_name) for sect in search_sect_obj ]

                        if new_const_value == None or len(new_const_value) == 0:
                            if allow_missing:
                                new_const_value = ""
                            else:
                                raise ValueError('Could not find keyword: %s in section: %s in file: %s' % (search_key_name, search_sect_name, extract_filename))
                        elif len(new_const_value) == 1:
                            new_const_value = new_const_value[0]

                        loaded_constants.append(new_const_name)
                        existingDict[ new_const_name ] = new_const_value

    if len(loaded_constants) > 0:
        logger.debug('Loaded values: %s' % ', '.join(loaded_constants))
                    
    return existingDict


def Read_Id_List_File(id_list_file, section=None, valuesDict={}, mapDict={}):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)
    
    if id_list_file == None:
        raise ValueError('id list file not defined')

    if section != None:
        logger.debug('Reading id list from section %s file: %s' % (section, id_list_file))
        id_list_str = Get_List_File_Values(id_list_file, str(id_list_file), section)
    else:
        logger.debug('Reading id list from file: %s' % id_list_file)
        id_obj = L2_Input.Input_File(id_list_file)                                   
        id_list_str = id_obj.Get_Matrix_Data()

    if id_list_str == None:
        #raise IOError('Did not read any sounding ids from id list file: %s' % id_list_file)
        return []

    id_list_long = []
    for curr_id_str in id_list_str:
        # Try to match sounding id pattern
        id_match = re.search('\d{3,17}\w?', curr_id_str)
        if not id_match:
            raise IOError('Could not find sounding id in string: "%s" in file %s' % (curr_id_str, id_list_file))
        beg_pos = int(id_match.start())
        end_pos = int(id_match.end())
        found_id = curr_id_str[beg_pos:end_pos]
        id_list_long.append( found_id )
        
    return id_list_long
    
def Parse_KeyVal_Str_List(str_list, orig_hash={}):

    if str_list != None:
        for sub_pair in str_list:
            if sub_pair == None:
                continue
            
            split_out = [ split_var.strip() for split_var in sub_pair.split('=', 1) ]
            if len(split_out) == 2:
                (sub_key, sub_value) = split_out
            else:
                (sub_key, sub_value) = (split_out, '')

            try:
                orig_hash[ sub_key ] = sub_value
            except TypeError:
                raise Exception('Error parsing key value string item: %s from list: %s' % (sub_pair, str_list))

    return orig_hash

def Set_Path_Constants(valueDict, controlFile):
    cf_path_dir   = os.path.realpath(os.path.dirname(os.path.expanduser(controlFile)))
    cf_path_parts = os.path.split(cf_path_dir)
    
    # Deprecate these eventually with better names...
    valueDict[ 'CONTROL_FILENAME' ]        = controlFile
    valueDict[ 'CONTROL_FILE_DIRECTORY' ]  = cf_path_dir
    valueDict[ 'CONTROL_FILE_PARENT_DIR' ] = cf_path_parts[0]

    return controlFile

def Process_Includes(sectionObj, valuesDict, allow_unresolved=False, recurseLevel=0):

    if recurseLevel > 100:
        raise RuntimeError('Indefinite recursive loop detected in includes')

    newChildren = []
    for childObj in sectionObj.children:
            
        if childObj.type == 'assignment' and childObj.leaf.lower() == 'include':
            if len(childObj.children) != 1:
                raise ValueError('Expected 1 value for include keyword, got %d' % len(childObj.children))
            
            include_tmpl = childObj.children[0].leaf

            try:
                include_file = Expand_Filename( Apply_Template(include_tmpl, valuesDict) )[0]
            except TemplateError:
                if allow_unresolved:
                    # Add back original include and continue onward
                    newChildren.append( childObj )
                    continue
                else:
                    raise

            if not os.path.exists(include_file):
                raise IOError('Included file "%s" does not exist' % include_file)

            includeObj = L2_Input.Input_File(include_file)

            # Change directory path constants to match that if included file
            inclValuesDict = copy.copy(valuesDict)
            Set_Path_Constants(inclValuesDict, include_file)

            Process_Includes(includeObj, inclValuesDict, allow_unresolved, recurseLevel+1)
            
            if includeObj.rootNode == None:
                raise IOError('Was unable to parse included file "%s"' % include_file)

            for include_child in includeObj.rootNode.children:
                newChildren.append(include_child)
            
        elif childObj.type == 'section':
            Process_Includes(childObj, valuesDict, allow_unresolved, recurseLevel+1)
            newChildren.append(childObj)
        else:
            newChildren.append(childObj)

    # Make sure we assign the children to the right place and don't get bitten
    # by the __getattr_ in L2_Input
    if hasattr(sectionObj, 'rootNode'):
        sectionObj.rootNode.children = newChildren
    elif hasattr(sectionObj, 'children'):
        sectionObj.children = newChildren
    else:
        raise Exception('Can not assign new children to section in object: %s' % sectionObj)

def HDF_Has_Dataset(hdf_filename, dataset_name):
    has_dataset = False
    with h5py.File(hdf_filename, 'r') as hdf_obj:
        try:
            dataset_values = hdf_obj[dataset_name]
            has_dataset = True
        except:
            pass

    return has_dataset
