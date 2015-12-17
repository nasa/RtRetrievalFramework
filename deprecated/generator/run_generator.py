#!/usr/bin/env python

from types import ListType
import os
import re
import sys
import copy
import shutil
import logging
from optparse import OptionParser

import L2_Input
from OCO_TextUtils import index_range_list, evaluate_bool_str
import L2_Log_Util

from Generator_Utils import *
from DiskOperations import Process_Operation_Section, COPYTREE_IGNORE_PATTERNS
import ThreadPool

LOADED_MODULES = {}

if os.path.exists('%s/%s' % (os.path.dirname(sys.argv[0]), 'modules')):
    sys.path.append('%s/%s' % (os.path.dirname(sys.argv[0]), 'modules'))

def Set_Combinations(value_lists_iter, level=0):
    """Creates a list of dictionaries for each combination of the items in a dictionary of lists
       For efficiency an iter object is required.
    """
    
    try:
        list_name, list_values = value_lists_iter.next()

        combination_list = []
        # For each sub combinations add our current list's values to that combination
        for sub_combination in Set_Combinations(value_lists_iter, level+1):
            if list_values != None:
                for curr_value in list_values:
                    new_combination = dict({list_name:curr_value}, **sub_combination)
                    combination_list.append( new_combination )
                
        return combination_list
    except StopIteration:
        # Empty combination
        return [ {} ]

def Process_Modules(parentSect, ignoreNameList, source_file, dest_file, valuesDict, fileMap, buffer_objs=None):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)
    
    moduleNames = []
    for parentSectName in parentSect.Get_All_Section_Names():
        sectNameParts = parentSectName.split('->')
        if len(sectNameParts) == 2:
            (baseSectName, moduleSectName) = sectNameParts

            if not moduleSectName in ignoreNameList and not moduleSectName in moduleNames:
                moduleNames.append(moduleSectName)

    if len(moduleNames) == 0:
        return

    # Store extra keywords in hash for possible module use
    parentKeywords = {}
    for key_name in parentSect.Get_All_Keyword_Names():
        key_val = Apply_Template( parentSect.Get_Keyword_Value(key_name), valuesDict, mapDict=fileMap )
        parentKeywords[key_name] = key_val
               
    for moduleName in moduleNames:

        moduleObj = None
        if moduleName in LOADED_MODULES:
            moduleObj = LOADED_MODULES[moduleName]
        else:
            tryModuleName = [moduleName.lower(), moduleName, moduleName.upper()]
            nameIdx = 0
            logger.debug('Loading module: %s' % moduleName)
            while moduleObj == None and nameIdx < len(tryModuleName):
                try:
                   moduleObj = __import__(tryModuleName[nameIdx], globals(), {}, [])
                except ImportError:
                    if str(sys.exc_value) == ('No module named %s' % tryModuleName[nameIdx]):
                        nameIdx += 1
                    else:
                        raise

            LOADED_MODULES[moduleName] = moduleObj
                
        if moduleObj == None:
            raise ImportError('Could not load module %s' % moduleName)

        allModuleSections = parentSect.Get_Section('->%s' % moduleName)

        # Use keyword that allows module to be rejected or accepted based on template evaluation
        matchedModuleSections = []
        for module_sect in allModuleSections:
            if Should_Process(module_sect, valuesDict, fileMap):
                matchedModuleSections.append(module_sect)

        # No matching modules, continue to the next modilename
        if len(matchedModuleSections) == 0:
            continue

        if buffer_objs != None:
            # Get a possibly pre-existing source object
            source_obj = buffer_objs.get(source_file, None)
            if source_obj == None and source_file != None:
                source_obj = NamedStringIO(filename=source_file)

            if source_obj != None:
                source_obj.seek(0)

            # Always write into a new object
            dest_obj = NamedStringIO()
            dest_obj.filename = dest_file
            dest_obj.seek(0)
        else:
            source_obj = source_file
            dest_obj   = dest_file

        try:
            logger.debug('Running module %s: %s -> %s' % (moduleName, source_obj, dest_obj))
            moduleArgs = (source_obj, dest_obj, parentKeywords, matchedModuleSections, valuesDict, fileMap, buffer_objs)
            moduleReturn = moduleObj.Process_File(*moduleArgs[:moduleObj.Process_File.func_code.co_argcount])

        except TemplateError:
            # Make template error more descriptive
            if (str(source_obj) != str(dest_obj)):
                logger.error('\nERROR: %s for source: %s, destination %s' % (sys.exc_value, source_obj, dest_obj))
            else:
                logger.error('\nERROR: %s for filename: %s' % (sys.exc_value, source))
            if hasattr(parentSect, 'filename'):
                logger.error('While processing script: %s' % parentSect.filename)
            return

        # if source == dest then this will overwrite a pre-existing object
        # but only do this if the output object was actually written to
        # ie. only if the file posistion is greater than the beginning
        # the modules should not seek to the beginning because they do not know
        # if they are appending to a larger file or not
        if buffer_objs != None and dest_obj.tell() > 0:
            buffer_objs[dest_file] = dest_obj

def Process_File_Sections(fileSectList, caseDir, valuesDict=None, mapDict=None):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)

    # Buffer file operations over all File Sections for a case directory
    buffer_objs = {}
    
    for fileSect in fileSectList:
        # The file map section overwrites the set map section if
        # it exits
        fileMapSections = fileSect.Get_Section('FILE->MAP')

        if len(fileMapSections) > 0:
            fileMap = Get_Map_Values( fileMapSections, valuesDict, existing=mapDict )
        else:
            # Use set map if not a file map defined
            fileMap = mapDict

        # Get filename to open and modify, make sure this is defined
        srcFilenameTmpl = fileSect.Get_Keyword_Value('source_filename')
        dstFilenameTmpl = fileSect.Get_Keyword_Value('dest_filename')

        if srcFilenameTmpl != None:
            srcModFilename = Expand_Filename( Apply_Template(srcFilenameTmpl, valuesDict, mapDict=fileMap) )
        else:
            srcModFilename = (None,)

        # If there is no different destination name then use source name(s)
        if dstFilenameTmpl == None or len(dstFilenameTmpl) == 0:
            dstFilenameVal = srcModFilename
        else:
            dstFilenameVal = Expand_Filename( Apply_Template(dstFilenameTmpl, valuesDict, mapDict=fileMap ) )

        if not hasattr(srcModFilename, '__iter__') or not hasattr(dstFilenameVal, '__iter__'):
            raise Exception('source filename: %s and destination %s must both be iterable' % (srcModFilename, dstFilenameVal))
            
        if len(srcModFilename) != len(dstFilenameVal):
            raise Exception('source filename: %s and destination %s must both the same length since they are both iterable' % (srcModFilename, dstFilenameVal))
            
        # Go through matching pairs of filenames performing the same operations
        for (curr_source, curr_dest) in zip(srcModFilename, dstFilenameVal):
            # Undefine the source file if it ends up being the same as the destination
            # and the source does not exist. Good chance it will be created into the destination
            # Otherwise an error would be thrown by a module expecting a file to be present
            if curr_source == curr_dest and curr_source != None and not os.path.exists(curr_source):
                curr_source = None

            internalSectNames = ['MAP']
            Process_Modules(fileSect, internalSectNames, curr_source, curr_dest, valuesDict, fileMap, buffer_objs)

    # If we were buffering objects then write all files that need writing
    # only destination filenames written into buffer dictionary so all of
    # them should be written
    if buffer_objs != None:
        for out_filename, buffer_obj in buffer_objs.items():
            buffer_obj.write_file(out_filename)
            buffer_obj.close_buffer()

def Process_Module_Dirs(module_dir_tmpls, valuesDict=None, mapDict=None):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)
    
    if not type(module_dir_tmpls) is ListType:
        module_dir_tmpls = [ module_dir_tmpls ]

    for curr_dir_tmpl in module_dir_tmpls:
        curr_dir = Apply_Template(curr_dir_tmpl, valuesDict, mapDict)

        if curr_dir != None and len(curr_dir) > 0 and not curr_dir in sys.path:
            curr_dir = Expand_Filename(curr_dir)[0]
            logger.debug('Adding modules dir: %s' % curr_dir)
            sys.path.append(curr_dir)

#######################
# Start of script code

def Process_Case(caseValDict, constDict, currControlSet, setMap, baseDir, removeExisting, ignored_module_names):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)
    
    # Print banner with list of value for case
    logger.info('-' * 80)
    if len(caseValDict.keys()) > 0:
        logger.info('Case dictionary contents:')
        maxKeywordLen = max([ len(x) for x in caseValDict.keys()])
        fmt_code = '%' + ('%d' % maxKeywordLen) + 's = %s'

        for (case_key, case_val) in caseValDict.iteritems():
            logger.info(fmt_code % (case_key, case_val))

        logger.info('')

    # Add non-preexisting constant substitutions to case substitutions dict
    for (sub_key, sub_val) in constDict.iteritems():
        if sub_key not in caseValDict:
            caseValDict[ sub_key ] = sub_val

    # Per case keywords. Use to over-ride global constants
    # Evaluate here after filling in path variables for the current case
    caseKeywordSects = currControlSet.Get_Section('->KEYWORDS')
    caseValDict = Get_Constant_Values(caseKeywordSects, existingDict=caseValDict, mapDict=setMap)

    # Create output dir from template
    sourceDir = Apply_Template( currControlSet.Get_Keyword_Value('source_dir'), caseValDict, mapDict=setMap )
    caseValDict[ 'source_dir' ] = sourceDir

    sourceFile = Apply_Template( currControlSet.Get_Keyword_Value('source_filename'), caseValDict, mapDict=setMap )
    caseValDict[ 'source_file' ] = sourceFile

    if currControlSet.Get_Keyword_Value('output_directory') != None:
        raise DeprecationWarning('output_directory keyword is deprecated')

    caseDir = Apply_Template( currControlSet.Get_Keyword_Value('dest_case_dir'), caseValDict, mapDict=setMap )
    caseValDict[ 'dest_case_dir' ] = caseDir

    # Any script directories per case
    if currControlSet.Get_Keyword_Value('scripts_dir') != None:
        raise DeprecationWarning('scripts_dir keyword is deprecated')

    Process_Module_Dirs(currControlSet.Get_Keyword_Value('modules_dir'), caseValDict, mapDict=setMap)

    skipIfExists = Apply_Template( currControlSet.Get_Keyword_Value('skip_if_exists'), caseValDict, mapDict=setMap )
    onlyIfExists = Apply_Template( currControlSet.Get_Keyword_Value('only_if_exists'), caseValDict, mapDict=setMap )

    # Expand paths should they contain ~/ or ~user
    sourceDir = Expand_Filename(sourceDir)[0]
    baseDir   = Expand_Filename(baseDir)[0]

    # Make sure required parameters
    if baseDir == None or len(baseDir) == 0:
        raise IOError('baseDir keyword must be defined for SET')

    # Create case dir filename from dest dir and output dir template
    if caseDir != None and len(caseDir) > 0:
        fullCaseDir = os.path.join(baseDir, caseDir)
        fullCaseDir = fullCaseDir.replace('//', '/') # Remove double // for astetics
    else:
        fullCaseDir = baseDir

    logger.info('Case directory:')
    logger.info(os.path.realpath(fullCaseDir))

    if skipIfExists != None and len(skipIfExists) > 0:
        if type(skipIfExists) is not ListType:
            skipIfExists = [ skipIfExists ]
        skipExistFiles = [ Expand_Filename(checkName, fullCaseDir)[0] for checkName in skipIfExists ]

        skipReason = None
        for currSkipCheck in skipExistFiles:
            if os.path.exists(currSkipCheck):
                skipReason = currSkipCheck
                break

        if skipReason != None:
            logger.debug('Skipping recreation of: %s because this file exist: %s' % (fullCaseDir, skipReason))
            return

    if onlyIfExists != None and len(onlyIfExists) > 0:
        if type(onlyIfExists) is not ListType:
            onlyIfExists = [ onlyIfExists ]
        onlyExistFiles = [ Expand_Filename(checkName, fullCaseDir)[0] for checkName in onlyIfExists ]
        skipReason = None
        for currExistCheck in onlyExistFiles:
            if not os.path.exists(currExistCheck):
                skipReason = currExistCheck
                break
        if skipReason != None:
            logger.debug('Skipping creation of: %s because this file does not exist: %s' % (fullCaseDir, skipReason))
            return

    if not Should_Process(currControlSet, caseValDict, setMap):
        logger.debug('Skipping creation of: %s because of template evaluation result' % (fullCaseDir))
        return

    # Make base part of case path if it does not exist
    # else remove the previous case directory
    if not os.path.exists(fullCaseDir):
        tcDirBase = (os.path.split( fullCaseDir.rstrip('/') ))[0]

        if len(tcDirBase) > 0 and not os.path.exists(tcDirBase):
            logger.debug('Creating case directory base: %s' % tcDirBase)
            os.makedirs(tcDirBase)

    elif sourceDir != None and len(sourceDir) > 0:
        # Only remove case directory if there is a source dir
        # defined. Would only need to remove for the purpose
        # of copying sourceDir
        if removeExisting:
            logger.debug('Removing old case directory: %s' % fullCaseDir)
            shutil.rmtree(fullCaseDir)
        else:
            raise IOError('Will not remove existing case directory: %s' % fullCaseDir)

    # If a source file is defined then copy its contents to the
    # case output dir
    if sourceDir != None and len(sourceDir) > 0:
        if not os.path.exists(sourceDir):
            errMsg =  'SET source_dir does not exists: %s' % sourceDir
            raise IOError(errMsg)

        logger.debug('Copying source data from %s to %s' % (sourceDir, fullCaseDir))
        try:
            shutil.copytree(sourceDir, fullCaseDir, symlinks=True, ignore=shutil.ignore_patterns(*COPYTREE_IGNORE_PATTERNS))
        except:
            raise IOError('Error copying from %s to %s' % (sourceDir, fullCaseDir))
    elif sourceFile != None and len(sourceFile) > 0:
        if not os.path.exists(sourceFile):
            errMsg =  'SET source_filename does not exists: %s' % sourceFile
            raise IOError(errMsg)
        try:
            shutil.copy(sourceFile, fullCaseDir)
        except:
            raise IOError('Error copying from %s to %s' % (sourceFile, fullCaseDir))
    else:
        # Since we did not copy a source directory, must create
        # case dir manually
        if not os.path.exists(fullCaseDir):
            logger.debug('Creating case directory: %s' % fullCaseDir)
            os.makedirs(fullCaseDir)

    # Change to case directory so that any scripts can reference relative files from here
    old_dir = os.getcwd()
    os.chdir(fullCaseDir)

    # Process setup for the current case. Performs actions such as
    # Deleting extraneous files and making necessary directories
    if len(currControlSet.Get_Section('->SETUP')) > 0:
        raise DeprecationWarning('SETUP section deprecated, use CASE_SETUP instead')

    for setupSection in currControlSet.Get_Section('->CASE_SETUP'):
        logger.debug('Executing CASE_SETUP section')

        Process_Operation_Section(setupSection, old_dir, caseValDict, setMap)

    # Loop over each file to be modified for the current list
    # values combination
    caseFileSects = currControlSet.Get_Section('->FILE')
    Process_File_Sections(caseFileSects, fullCaseDir, caseValDict, setMap)

    internalSectNames = ['FILE', 'CASE_SETUP', 'SET_SETUP', 'CLEANUP', 'LIST', 'MAP', 'KEYWORDS']

    if ignored_module_names != None:
        internalSectNames += ignored_module_names

    Process_Modules(currControlSet, internalSectNames, None, fullCaseDir, caseValDict, setMap)

    # Process cleanup for the current case. Performs actions such as
    # Deleting extraneous files and making necessary directories
    for cleanupSection in currControlSet.Get_Section('->CLEANUP'):
        logger.debug('Executing CLEANUP section')

        Process_Operation_Section(cleanupSection, old_dir, caseValDict, setMap)

    # Change back to previous directory
    os.chdir(old_dir)
   

def run_generator(controlFileList, constDictIn={}, maxNumberCases=None, removeExisting=False, export_path=None, ignored_module_names=None):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)
   
    # Run generator for each specified file
    for fileIdx, controlFile in enumerate(controlFileList):

        # Reset for every control file
        constDict = copy.copy(constDictIn)

        logger.info('=' * 80)
        logger.info('Processing: %s' % controlFile)
        logger.info('')

        # Add helper substititions
        Set_Path_Constants(constDict, controlFile)

        # Load control file
        if not os.path.exists(controlFile):
            logger.error('Control file: %s does not exist' % controlFile)
            return

        controlObj = L2_Input.Input_File(controlFile)

        # Make sure file was parsed
        if controlObj.rootNode == None:
            raise IOError('Was unable to parse %s' % controlFile)

        # Load include files that may include CONSTANT blocks
        # but do not die if there are unresolved template replacements
        Process_Includes(controlObj, constDict, allow_unresolved=True)

        # Load constants blocks into constants dictionary
        constSections = controlObj.Get_Section('CONSTANTS')
        constDict = Get_Constant_Values(constSections, existingDict=constDict)

        # Load remaining includes but make sure there are no unresolved replacements
        Process_Includes(controlObj, constDict, allow_unresolved=False)

        # Process file wide modules
        Process_Module_Dirs(controlObj.Get_Keyword_Value('modules_dir'), constDict)

        # Export control file after all includes have been loaded
        if export_path != None:
            orig_file_comment = L2_Input.Node('comment', leaf='# Original file: %s' % os.path.realpath(controlFile))
            controlObj.rootNode.children.insert(0, orig_file_comment)
            controlObj.Write(os.path.join(export_path, '%d-%s' % (fileIdx, os.path.basename(controlFile))))

        # Look for old named TESTSET section
        if len(controlObj.Get_Section('TESTSET')) > 0:
            raise DeprecationWarning('TESTSET name for groups is deprecated, use name SET instead')

        # Find all sets in the control file
        controlSets = controlObj.Get_Section('SET')

        # Loop over each set and create cases based on description in file
        for currControlSet in controlSets:
            # Get set map values
            mapSections = currControlSet.Get_Section('SET->MAP')
            setMap = Get_Map_Values(mapSections, constDict)

            # Get list names and associate values
            listSections = currControlSet.Get_Section('SET->LIST')
            listValuesDict = Get_List_Section_Values(listSections, constDict, setMap)

            # Get set specific template values
            if currControlSet.Get_Keyword_Value('dest_dir') != None:
                raise DeprecationWarning('dest_dir keyword is deprecated')

            baseDir = Apply_Template( currControlSet.Get_Keyword_Value('dest_base_dir'), constDict, mapDict=setMap )

            # Add useful items to constDict 
            constDict[ 'dest_base_dir' ] = baseDir
            constDict[ 'overwrite' ] = str(removeExisting)

            for setupSection in currControlSet.Get_Section('->SET_SETUP'):
                logger.debug('Executing SET_SETUP section')

                # Merge contants dictionary with dictionary of list values so that
                # all possibilities for a given LIST can be matched for file
                # operations can be done before starting case operations                
                setSetupDict = dict(constDict, **listValuesDict)
                Process_Operation_Section(setupSection, os.curdir, setSetupDict, setMap)

            # Create combinations matrix given list of list values
            caseCombinations = Set_Combinations(listValuesDict.iteritems())

            ##################################################################
            # Loop over each combination of list values and create a case
            # for each
            if maxNumberCases != None:
                combMatrix = combMatrix[:min(maxNumberCases+1, len(combMatrix))]

            caseArgs = []
            for caseValDict in caseCombinations:
                if not caseValDict in caseArgs:
                    caseArgs.append(caseValDict)
                else:
                    raise Exception('Duplicate caseValDict: %s' % caseValDict)

                Process_Case(caseValDict, constDict, currControlSet, setMap, baseDir, removeExisting, ignored_module_names)

            def Process_Case_Closure(caseValDict):
                # None of these closured arguments should be modified at all by called routine
                return Process_Case(caseValDict, constDict, currControlSet, setMap, baseDir, removeExisting, ignored_module_names)
                
            #mt = ThreadPool.MultiThread( Process_Case_Closure, caseArgs, maxThreads=1 )
            #mt.start()
            #mt.join()


##########################

if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog [options] <control_file> [ <control_file> ... ]")

    parser.add_option( "-m", "--max_cases", dest="max_cases",
                       metavar="NUM",
                       type="int",
                       help="maximum number of case to generate for each control file"
                       )

    parser.add_option( "-c", "--constant", dest="constants",
                       metavar="KEY=VALUE",
                       type="string",
                       action="append",
                       help="define a constant keyword substitution and value"
                       )

    parser.add_option( "-o", "--overwrite", dest="overwrite",
                       default=False,
                       action="store_true",
                       help="Overwrite existing case directories"
                       )

    # Parse command line arguments
    (cmd_options, control_files) = parser.parse_args()

    # Initialize logging
    L2_Log_Util.init_logging()

    if len(control_files) == 0:
        parser.error('at least one testcase input file must be specified')

    # Create list substitution dictionary from command line args
    constants_dict = Parse_KeyVal_Str_List(cmd_options.constants)

    if cmd_options.max_cases != None:
        max_num_cases = cmd_options.max_cases
    else:
        max_num_cases = None

    run_generator(control_files, constants_dict, max_num_cases, cmd_options.overwrite)
