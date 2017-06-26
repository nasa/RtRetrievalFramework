import os
import copy
import logging

from Generator_Utils import *
from OCO_TextUtils import evaluate_bool_str

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    logger = logging.getLogger(os.path.basename(__file__))    

    optionDict = copy.copy(valuesDict)
    for keyName in fileKeywords.keys():
        optionDict[keyName] = fileKeywords[keyName]

    if source != None:
        logger.warning('source ignored by EXECUTE module')

    if destination != None:
        logger.warning('destination ignored by EXECUTE module')

    for executeSect in moduleSections:
        binary  = Apply_Template(executeSect.Get_Keyword_Value('binary'), valuesDict, mapDict=mapDict)
        options = Apply_Template(executeSect.Get_Keyword_Value('options'), optionDict, mapDict=mapDict)
        chdir   = Apply_Template(executeSect.Get_Keyword_Value('chdir'), optionDict, mapDict=mapDict)
        quiet   = evaluate_bool_str(executeSect.Get_Keyword_Value('quiet'))

        if binary == None or len(binary) == 0:
            raise ValueError("No binary name specified")

        if type(options) is ListType:
            options = ' '.join(options)

        old_dir = os.getcwd()
        if chdir != None:
            if not os.path.exists(chdir):
                raise IOError('Could not change to dir: %s' % chdir)
            os.chdir(chdir)

        envSects = executeSect.Get_Section('EXECUTE->ENVIRONMENT')
        if envSects != None:
            logger.debug('Setting enviromental variables: ')
            for currEnvSect in envSects:
                for keyName in currEnvSect.Get_All_Keyword_Names():
                    os.environ[keyName] = Apply_Template(currEnvSect.Get_Keyword_Value(keyName), valuesDict, mapDict=mapDict)
                    logger.debug(keyName, '=', os.environ[keyName])

        if options != None:
            run_command = '%s %s' % (binary, options)
        else:
            run_command = '%s' % (binary)

        logger.debug('Executing command: %s' % run_command)

        run_obj = os.popen(run_command)
        if not quiet:
            for run_line in run_obj.readlines():
                logger.debug(run_line.strip())
        run_obj.close()

        os.chdir(old_dir)

