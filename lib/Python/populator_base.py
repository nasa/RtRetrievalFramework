from __future__ import division
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import os
import logging
import glob
import sys
import re
import six
from full_physics.l2_input import L2InputFile
from abc import ABCMeta, abstractmethod, abstractproperty

CLUSTER_OPTIONS = {
        'torque': {
            'array_idx_var_name': 'PBS_ARRAYID',
            'job_array_option': '-t',
            'depend_option': '-W depend=afteranyarray'
            },
        'pbs_pro': {
            'array_idx_var_name': 'PBS_ARRAY_INDEX',
            'job_array_option': '-J',
            'depend_option': '-W depend=afterany'
             },
        'pleiades': {
            'array_idx_var_name': 'PBS_ARRAY_INDEX',
            'job_array_option': '-J',
            'depend_option': '-W depend=afterany'
            },
        }

class PopulatorBase(object):
    '''This is the base class for various Populator classes
    (e.g., GosatPopulator). This contains the routines that are common
    to all the populators.

    A populator create a number of files from templates, filling in values
    based on defaults, calculations, and user supplied information. These
    files can then be used to actually run Level 2 Full Physics to match
    the selected values.
    '''

    def __init__(self, **user_settings):
        '''Initialize object'''
        self.logger = logging.getLogger(os.path.basename(__file__))
        self.aggregate = False
        self.skip_check = False
        self.l2_binary_filename = None
        self.abscoversion = ""
        self.group_size = 1
        self.use_subdirectory = False
        self.parallel_size = 1
        self.email_address = ""

        # Add cluster specific template substition values
        if 'target_cluster' in user_settings:
            cluster_settings = CLUSTER_OPTIONS.get(user_settings['target_cluster'], None)
            if cluster_settings == None:
                raise ValueError('Target cluster name: "%s" is unrecognized, use on of the following: %s' % (user_settings['target_cluster'], list(CLUSTER_OPTIONS.keys())))
            user_settings.update(cluster_settings)

        # Add dictionary of external settings and substitutions as 
        # attributes of the class
        for k, v in list(user_settings.items()):
            setattr(self, k, v)

        # A derived class can make this false in order to suppress trying
        # to create a l2/l1 aggregate. This is a useful thing for all the ACOS
        # like populators, but doesn't make any sense for FTS.
        self.have_l1b1 = True
        self._l2_input_file_cache = {}

    # This will get filled in by a map of type name and class to handle
    # it, e.g., populator_list["oco"] = OcoPopulator. We fill this in as each
    # class gets loaded, so we automatically know this without needing to
    # create a central list.
    populator_list = { }

    @abstractproperty
    def config_file_obj_re(self):
        '''Regular expression used to identify a configuration file object
        as a particular type (e.g., FTSFullPhysics for FtsPopulator) '''
        pass

    @abstractproperty
    def config_template_file(self):
        '''This should return the template file that create_config.py should use
        to create the configuration file.'''
        pass

    @abstractproperty
    def variable_exports_list(self):
        '''List of variables and values to include in run file.'''
        pass

    @abstractproperty
    def additional_var_list(self):
        '''List of variables to include in fun file. This is indexed by
        job index, unlike the variable_exports_list'''
        pass

    @abstractmethod
    def populate(self, config_filename, **user_settings):
        '''This is the function that takes the configuration file and
        user settings, and use this to generate the run scripts.'''
        pass

    @property
    def input_config_filename(self):
        '''Configuration file. Note that self.processing_dir should be
        filled in first before calling this.'''
        return self.processing_dir + "/sdos_input_list.dat"

    @property
    def sounding_id_list_filename(self):
        return self.processing_dir + "/sounding_id.list"
    @property
    def job_script_filename(self):
        return self.processing_dir + "/l2_fp_job.sh"
    @property
    def pleiades_job_script_filename(self):
        return self.processing_dir + "/pleiades_job.sh"
    @property
    def aggregate_script_filename(self):
        return self.processing_dir + "/aggregate.sh"
    @property
    def launch_script_filename(self):
        return self.processing_dir + "/launch_jobs.sh"
    @property
    def log_directory(self):
        return self.processing_dir + "/log"

    # The next set of properties are all things needed by the templates. They
    # just request the value when it is needed, and this class calculates it.
    @property
    def output_directory(self):
        return self.processing_dir + "/output"
    @property
    def python_path(self):
        return os.environ["PYTHONPATH"]
    @property
    def abscodir(self):
        return os.environ.get("abscodir", "/groups/algorithm/l2_fp/absco")
    @property
    def merradir(self):
        return os.environ.get("merradir", "/groups/algorithm/l2_fp/merra_composite")
    @property
    def bin_path(self):
        return os.environ.get("PATH", "")
    @property
    def lib_path(self):
        return os.environ.get("LD_LIBRARY_PATH", "")
    @property
    def l2_lua_path(self):
        res = os.environ["LUA_PATH"]
        if self.l2_config_filename:
            user_config_dir = os.path.realpath(os.path.dirname(self.l2_config_filename))
            res = user_config_dir + "/?.lua;" + res
        return res
    @property
    def l2_support_path(self):
        return os.environ.get("L2_SUPPORT_PATH", "")
    @property
    def l2_support_util_path(self):
        return os.environ.get("L2_SUPPORT_UTIL_PATH", "")
    @property
    def qsub_name(self):
        return "l2_fp"
    @property
    def aggregate_name(self):
        return "aggregate"
    @property 
    def addl_var_list_def(self):
        return "\n".join(
            ['%s_list=( `cat %s` )\n' % (k, v) 
             for k, v in list(sorted(self.additional_var_list.items()))])
    @property
    def addl_var_run_export(self):
        return "\n".join(
            ['export %s="${%s_list[${job_index}]}"\n' % (i, i)
             for i in list(sorted(self.additional_var_list.keys()))])
    @property
    def variable_exports(self):
        return "\n".join(
            ['export %s="%s"\n' % (var_name, var_value) 
             for var_name, var_value in list(sorted(self.variable_exports_list.items()))])
    @property
    def spectrum_file(self):
        return self.variable_exports_list["spectrum_file"]

    @staticmethod
    def create_populator_from_config_type(config_type, **user_settings):
        '''This searches the populator_list for the given config type, and
        if found returns an populator object for that config type. If we
        don't find anything, we return None.'''
        if config_type in PopulatorBase.populator_list:
            return PopulatorBase.populator_list[config_type](**user_settings)
        else:
            return None

    def _l2_input_file(self, file):
        '''This caches reading a config file, so we don't parse the same file 
        multiple_times.'''
        if(file not in self._l2_input_file_cache):
            self._l2_input_file_cache[file] = L2InputFile(file)
        return self._l2_input_file_cache[file]
    
    @staticmethod
    def create_populator_from_config_file(config_file, **user_settings):
        '''Read the L2 input configuration file supplied, and based on the
        type return the populator object that matches the type. If we can't
        find it, return None.'''
        l2_config = L2InputFile(config_file)
        t = l2_config.get_section("input->PGENameGroup")
        if len(t) > 0:
            pge = t[0].get_keyword_value("PGEName")
            if(pge is not None and len(pge) > 0):
                for k, v in list(PopulatorBase.populator_list.items()):
                    t = v(**user_settings)
                    if(re.search(pge, t.config_file_obj_re)):
                        return t
        
        # Some older config files don't have the PGE in it, so try working
        # off of the file name
        for k, v in list(PopulatorBase.populator_list.items()):
            if(re.search(config_file, r'^%s_' % k)):
                return v(**user_settings)
        return None
        
    def create_run_scripts(self, sounding_id_file, sounding_id_sect):
        '''Using the set value and the templates, create run scripts'''

        id_list = self.read_id_list_file(sounding_id_file, sounding_id_sect)
    
        if len(id_list) > self.group_size:
            self.max_array_index = old_div(len(id_list), self.group_size)
        else:
            # There is really only 1 index to be run, but without specficially setting this to 1
            # the range would be 0-0. But qsub does not like that range and will complain. So
            # we give it one more index so it at least runs the one we want and run_job will just 
            # skip the extra index
            #
            # Edge condition of either 1 sounding id or a single group of soundings smaller than
            # the grouping size
            self.max_array_index = 1

        # Subtract 1 if the group size is an integer multiple of the number of ids
        if self.max_array_index * self.group_size == len(id_list):
            self.max_array_index -= 1
        
        with open(self.sounding_id_list_filename, 'w') as si_fo:
            si_fo.write( ' '.join(id_list) )

        # Read in the templates. These are complicated enough we store them
        # in a separate file.
        tmpl_dir = os.path.dirname(__file__) + "/template/"
        if(self.have_l1b1):
            agg_template= open(tmpl_dir + "aggregate.tmpl").read()
        else:
            agg_template= open(tmpl_dir + "aggregate_nol1b1.tmpl").read()
        
        with open(self.aggregate_script_filename, 'w') as agg_fo:
            agg_fo.write( agg_template.format(self) )
        os.chmod(self.aggregate_script_filename, 0o755)
        
        job_template = open(tmpl_dir + "job.tmpl").read()
        with open(self.job_script_filename, 'w') as job_fo:
            job_fo.write( job_template.format(self) )
        os.chmod(self.job_script_filename, 0o755)
            
        if self.target_cluster == 'pleiades':
            pleiades_template = open(tmpl_dir + "pleiades_job.tmpl").read()
            with open(self.pleiades_job_script_filename, 'w') as pjob_fo:
                pjob_fo.write( pleiades_template.format(self) )
            os.chmod(self.pleiades_job_script_filename, 0o755)
        else:
            launch_template = open(tmpl_dir + "launch_jobs.tmpl").read()
            empty_launch_template = open(tmpl_dir + "empty_launch_jobs.tmpl").read()
            with open(self.launch_script_filename, 'w') as launch_fo:
                if len(id_list) == 0:
                    launch_fo.write( empty_launch_template.format(self) )
                else:
                    launch_fo.write( launch_template.format(self) )
            os.chmod(self.launch_script_filename, 0o755)

    def __expand_filename(self, filename, basePath=None):
        '''Expand a filename, used internally. This expands ~, and any
        wildcards'''
        if filename == None:
            return [ '' ]
    
        # Expand any usage of ~/ or ~user in pathnames
        filename = os.path.expanduser(filename)
        if basePath != None:
            basePath = os.path.expanduser(basePath)
   
        if(basePath != None and not filename[0] == '/'
           and not filename.find('..') == 0):
            fullname = basePath + '/' + filename
        else:
            fullname = filename

        matches = glob.glob(fullname)

        if len(matches) > 0:
            return matches
        else:
            return [ fullname ]

    def __get_list_file_values(self, listLocation, listName, sectionName=None, directoryLevels=None):
        # Expand any globs
        if isinstance(listLocation, six.string_types):
            listLocation = self.__expand_filename(listLocation)
    
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
           
            if isinstance(listFile, six.string_types) and os.path.isdir(listFile):
                self.logger.debug('Loading LIST %s contents from directory: %s' % (listName, listFile))
    
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
                self.logger.debug('Loading LIST %s contents from file: %s' % (listName, listFile))

                if isinstance(listFile, six.string_types):
                    listFileObj = open(listFile, 'r')
                elif hasattr(listFile, 'read'):
                    listFileObj = listFile
                else:
                    raise Exception('Unknown read object: %s' % listFile)

                for listLine in listFileObj.readlines():
                    if len(listLine.strip()) > 0 and listLine.strip()[0] != '#':
                        fileValues.append(listLine.strip())

                if isinstance(listFile, six.string_types):
                    listFileObj.close()
            else:
                self.logger.debug('Loading LIST %s section as %s contents from file: %s' % (sectionName, listName, listFile))
                fileObj = self._l2_input_file(listFile)
                sectNameParts = sectionName.split('->')

                foundSects = fileObj.get_section('->'.join(sectNameParts[0:-1]) + '->LIST')

                fileListSect = None
                for currFileSect in foundSects:
                    currListName = currFileSect.get_keyword_value('name')
                    if currListName != None and currListName == sectNameParts[-1]:
                        fileListSect = currFileSect.get_section('LIST->VALUES')
                        break
                if fileListSect == None or len(fileListSect) == 0:
                    raise IOError('Could not find section %s in file: %s' % (sectionName, listFile))

                fileValues = fileListSect[0].get_matrix_data()

        return fileValues

    def read_id_list_file(self, id_list_file, section=None):
        '''Read the sounding ids from the given list of files. The file
        list can contain wildcards and/or ~'''

    # Use old L2_Input syntax
        if section != None:
            section = section.replace("/", "->")

        if section != None:
            self.logger.debug('Reading id list from section %s file: %s' % (section, id_list_file))
            id_list = self.__get_list_file_values(id_list_file, str(id_list_file), section)
        else:
            self.logger.debug('Reading id list from file: %s' % id_list_file)
            # Quicker read for text only file
            id_list = open(id_list_file).read().split()

        if id_list == None:
            return []

        # Remove any white space
        id_list = [i.strip() for i in id_list]
        # Check for any bad data
        bad = [i for i in id_list if not re.match('\d{3,17}', i)]
        if(len(bad) > 0):
            raise IOError('Bad data in sounding ID list, could not find sounding id in string: "%s" in file %s' % (bad[0], id_list_file))
        return id_list

    def get_config_keyword_value(self, config_filename, keyword_path):
        '''Read a L2 input file as keyword/value pairs, and return the value
        for the given keyword'''
        config_obj = self._l2_input_file(config_filename)

        search_sect_name = '->'.join(keyword_path.split('/')[0:-1])
        search_key_name  = keyword_path.split('/')[-1]

        search_sect_obj = config_obj.get_section(search_sect_name)

        if search_sect_obj == None or len(search_sect_obj) == 0:
            raise KeyError('Could not find section: %s in file: %s' % (search_sect_name, config_filename))

        found_values = [ sect.get_keyword_value(search_key_name) for sect in search_sect_obj ]

        if len(found_values) == 1:
            return found_values[0]
        else:
            return found_values


    def is_processable(self, sounding_id_file, sounding_id_sect):
        '''Determine if we have any sounding ids to process, and if
        not log a warning message. This is also a hook if we have some
        more complicated test in the future (e.g., check that all the input
        data is available and readable. This return True if data can be
        processed, False otherwise (currently this is always True)'''
        # Load sounding id list to see if there are any soundings to process
        sounding_id_list = self.read_id_list_file(sounding_id_file, sounding_id_sect)

        if len(sounding_id_list) == 0:
            # Warning message, but we still want to process. It obviously doesn't
            # make a lot of sense to process 0 soundings, but handling this
            # ends up being a useful edge case for SDOS.
            self.logger.info("%s has no sounding ids, but proceeding anyways" % (sounding_id_file))
            return True
        else:
            self.logger.debug("%s is processable and has %d sounding ids" % (sounding_id_file, len(sounding_id_list)))
            return True

    def init_processing_dir(self):
        '''Initialize the processing directory, creating the various output
        directories'''
        if os.path.exists(self.input_config_filename):
            self.logger.debug("Removing input_config_filename: %s" % self.input_config_filename)
            os.remove(self.input_config_filename)
    
        create_dirs = ( self.log_directory,
                        self.output_directory, )
        for new_dir in create_dirs:
            if not os.path.exists(new_dir):
                self.logger.debug("Creating directory: %s" % new_dir)
                os.makedirs(new_dir)

# Miscellenous function used by both populate and create_config. We might
# remove this later, but for now leave here.

def parse_keyval_str_list(str_list, orig_hash={}):

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
