import os
import sys
import shutil
import stat
import logging

from OCO_TextUtils import index_range_list, evaluate_bool_str
from Generator_Utils import *

# Do not copy these patterns since they are not useful in the target of the copies
COPYTREE_IGNORE_PATTERNS = ('.svn', '*~', '.git_create_me', '.gitignore')

def Perform_DELETE(source, options):
    if os.path.isdir(source) and not os.path.islink(source):
        shutil.rmtree(source)
    else:
        os.remove(source)

def Perform_MKDIR(destination_dir, options):
    os.makedirs(destination_dir)

def Perform_MOVE(source, destination, options):
    shutil.move(source, destination)

def Perform_COPY(source, destination, options):       
    if os.path.isdir(source) and not os.path.islink(source):
        shutil.copytree(source, destination, ignore=shutil.ignore_patterns(*COPYTREE_IGNORE_PATTERNS))
    else:
        shutil.copy(source, destination)
        file_stats = os.stat(destination)
        os.chmod(destination, file_stats.st_mode | stat.S_IWRITE)

def Perform_LINK(source, destination, options):
    link_relative = evaluate_bool_str( options.get('link_relative'), False )
    remove_broken_link = evaluate_bool_str( options.get('remove_broken_link'), False )
    
    if link_relative:
        source = os.path.relpath(source, os.path.dirname(destination))
        
    os.symlink(source, destination)

    # Check for broken link
    if not os.path.exists(destination):
        # Remove broken symbolic link if requested
        if remove_broken_link:
            os.remove(destination) 
        raise IOError('Could not create or broken symbolic link %s named %s' % (source, destination))
    
def Process_Operation_Section(fileOpSection, baseSourceDir=None, valuesDict=None, mapDict=None):
    logger = logging.getLogger(sys._getframe(0).f_code.co_name)
   
    operation_options = fileOpSection.Get_Keywords_Dict()
    for curr_key, curr_val in operation_options.items():
        operation_options[curr_key] = Apply_Template(curr_val, valuesDict, mapDict=mapDict)

    fail_on_error   = evaluate_bool_str( operation_options.get('fail_on_error'), True)
    skip_if_exists  = operation_options.get('skip_if_exists')
    remove_existing = evaluate_bool_str( operation_options.get('remove_existing'), False )
   
    if skip_if_exists != None and len(skip_if_exists) > 0:
        skipExistFiles = Expand_Filename(skip_if_exists)
        if os.path.exists(skipExistFiles[0]):
            logger.debug('Skipping %s section because this file exist: %s' % (fileOpSection.leaf[0], skipExistFiles[0]))
            return

    if not Should_Process(fileOpSection, valuesDict, mapDict):
        logger.debug('Skipping %s section because of template evaluation result' % (fileOpSection.leaf[0]))
        return

    for fileAction in fileOpSection.Get_All_Section_Nodes():

        actionName = fileAction.leaf[0].upper()
        action_matrix_data = fileAction.Get_Matrix_Data()

        # Allow for a replacement map inside of file operation section
        # Trickily the order of the MAP can matter!
        if actionName == 'MAP':
            mapDict = Get_Map_Values( fileAction, existing=mapDict )

        elif fileOpSection.leaf[0] != actionName:
            try:
                actionFunc = eval('Perform_%s' % actionName)
            except NameError:
                raise NameError('Could not find function to handle %s operation %s' % (fileOpSection.leaf[0], actionName))

            for templateLine in action_matrix_data:                
                if hasattr(templateLine, '__iter__'):
                    sources      = Expand_Filename( Apply_Template( templateLine[0], valuesDict, mapDict=mapDict ), baseSourceDir )
                    destinations = Expand_Filename( Apply_Template( templateLine[1], valuesDict, mapDict=mapDict ) )
                else:
                    sources      = (None,)
                    destinations = Expand_Filename( Apply_Template( templateLine, valuesDict, mapDict=mapDict ) )
                
                for curr_source, curr_destination in zip(sources, destinations):
                   
                    try:
                        if os.path.exists(curr_destination) and remove_existing:
                            logger.debug('Removing existing destination: %s' % curr_destination)
                            Perform_DELETE(curr_destination, operation_options)

                        if curr_source == None:
                            logger.debug('Performing disk operation %s for section %s with file: %s' % (actionName, fileOpSection.leaf[0], curr_destination))
                            actionFunc(curr_destination, operation_options)
                        else:
                            logger.debug('Performing disk operation %s for section %s with source: %s, destination: %s' % (actionName, fileOpSection.leaf[0], curr_source, curr_destination))
                            actionFunc(curr_source, curr_destination, operation_options)
                            
                    except:
                        err_msg = 'Could not process disk operation: %s for section: %s with source: %s, destination: %s' % (actionName, fileOpSection.leaf[0], curr_source, curr_destination)
                        if fail_on_error:
                            logger.error(err_msg)
                            raise
                        else:
                            logger.debug(err_msg)




