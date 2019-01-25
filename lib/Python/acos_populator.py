from __future__ import absolute_import
from .populator_base import PopulatorBase
from full_physics.version_util import source_version, binary_version
import six
from full_physics.oco_matrix import OcoMatrix
import os

class AcosPopulator(PopulatorBase):
    '''This is a populator that has what is common for the Acos like
    populators.
    '''

    def __init__(self, config_sounding_id_section, config_input_keywords,
                 **user_settings):
        PopulatorBase.__init__(self, **user_settings)
        self.config_sounding_id_section = config_sounding_id_section
        self.config_input_keywords = config_input_keywords


    @property
    def variable_exports_list(self):
        return self.input_filenames

    @property
    def additional_var_list(self):
        return { }

    def set_input_config_values(self, sounding_id_file, sounding_id_sect, input_config_filename, input_file_list=[], **kwargs):
        '''This writes the input configuration values to the input configuration
        file (i.e., the sdos_input_list.dat file)'''

        keyword_defaults = {
            'file_id':                   'Scalar Retrieval Outputs',
            'exe_path':                  None,
            'exe_version':               None,
            'data_version':              None,
            'release_version':           None,
            'comments':                  '',
            'algorithm_descriptor':      None,
            'algorithm_maturity':        None,
            'l2_input_path':             None,
            'number_soundings':          None,
            }
        # Load sounding id list so we can leave a count in the config file we produce
        sounding_id_list = self.read_id_list_file(sounding_id_file, sounding_id_sect)

        file_keywords = {}
        file_keywords.update(keyword_defaults)
        file_keywords.update(kwargs)

        file_keywords['number_soundings'] = len(sounding_id_list)

        # Try getting versions from the binary file itself first
        exe_version = None
        data_version = None
        if 'exe_path' in file_keywords and file_keywords['exe_path'] != None:
            try:
                ver_ret = binary_version(file_keywords['exe_path'])
            except OSError as exc:
                raise OSError("Could not execute L2 binary: %s due to error: %s" % (file_keywords['exe_path'], exc))

            if ver_ret:
                file_keywords['release_version'] = ver_ret[0]
                exe_version = ver_ret[1]
                self.logger.debug('Retrieved release_version "%s" from binary %s' % (file_keywords['release_version'], file_keywords['exe_path']))
                self.logger.debug('Retrieved exe_version "%s" from binary %s' % (exe_version, file_keywords['exe_path']))

                data_version = ver_ret[2]
                if data_version != None:
                    self.logger.debug('Retrieved data_version "%s" from binary %s' % (data_version, file_keywords['exe_path']))

            # If the binary doesn't have any version information, then try looking for a CM directory
            # where the executable lives
            if exe_version == None:
                exe_dir = os.path.dirname(file_keywords['exe_path'])
                exe_version = source_version(exe_dir)
                if exe_version != None:
                    self.logger.debug('Retrieved exe_version "%s" from binary containing directory %s' % (exe_version, exe_dir))

        # If the binary is not in a source controlled directory try the src_path, which probably
        # came from an enviromental variable
        if exe_version == None and 'src_path' in file_keywords and file_keywords['src_path'] != None:
            exe_version = source_version(file_keywords['src_path'])
            if exe_version != None:
                self.logger.debug('Retrieved exe_version "%s" from source directory %s' % (exe_version, file_keywords['src_path']))

        if exe_version:
            file_keywords['exe_version'] = exe_version
        else:
            self.logger.error("Could not determine exe_version from executable: %s or source path: %s" % (file_keywords['exe_path'], file_keywords['src_path']))

        # If there was no binary version extracted from the binary then search for it from
        # the data_path variable
        if data_version != None:
            file_keywords['data_version'] = data_version
        elif 'data_path' in file_keywords and file_keywords['data_path'] != None:
            data_version = source_version(file_keywords['data_path'])
            self.logger.debug('Retrieved data_version "%s" from %s' % (file_keywords['data_version'], file_keywords['data_path']))

        if data_version:
            file_keywords['data_version'] = data_version
        else:
            self.logger.error("Could not determine data_version from path: %s" % file_keywords['data_path'])

        if 'L2_INPUT_PATH' in os.environ:
            file_keywords['l2_input_path'] = os.environ['L2_INPUT_PATH']

        self.logger.debug('Writing input file config file: %s' % input_config_filename)
        out_mat_obj = OcoMatrix()

        # Set items into input config file from values specified in configuraiton file
        for head_key_name, head_key_value in file_keywords.items():
            if hasattr(out_mat_obj, head_key_name):
                self.logger.debug('Set %s as an attribute' % head_key_name)

                prev_value = getattr(out_mat_obj, head_key_name)
                setattr(out_mat_obj, head_key_name, head_key_value)
            else:
                self.logger.debug('Set %s into header' % head_key_name)

                if isinstance(head_key_value, six.binary_type):
                    head_key_value = head_key_value.decode('UTF-8')

                if isinstance(head_key_value, six.string_types) and head_key_value.find(' ') >= 0:
                    out_mat_obj.header[head_key_name] = '"%s"' % head_key_value
                elif head_key_value == None:
                    out_mat_obj.header[head_key_name] = 'VALUE NOT SET'
                else:
                    out_mat_obj.header[head_key_name] = '%s' % head_key_value

        out_mat_obj.data = [fn for fn in input_file_list if fn != None and len(fn) > 0]
        out_mat_obj.write(input_config_filename, auto_size_cols=False)

    def populate(self, config_filename):
        '''This is the function that takes the configuration file and
        user settings, and use this to generate the run scripts.'''
        self.processing_dir = os.path.dirname(os.path.abspath(config_filename))

        # Check if the configuration can be processed, else signal
        # that we can not proceed
        if not self.skip_check:
            self.logger.info("Checking if processable")
            if not self.is_processable(config_filename,
                                       self.config_sounding_id_section):
                return False

        # Create common necessary files and directories
        self.logger.info("Initializing processing dir")

        # Get filenames for files used in processing from config file
        # Use Python 2.6 safe comprehension
        self.input_filenames = {}
        for inp_file_k,inp_file_v in list(self.config_input_keywords.items()):
            # If the section isn't in the input file (e.g., AlternativeSettings),
            # then treat this as just not finding the keyword.
            try:
                config_val = self.get_config_keyword_value(config_filename, inp_file_v)
            except KeyError:
                config_val = None
            if config_val != None and config_val == "NONE":
                config_val = None
            if config_val != None and len(config_val) > 0:
                self.input_filenames[inp_file_k] = os.path.realpath(config_val)
            elif inp_file_v.find("AlternativeSettings") < 0:
                # Only put empty declaration if file name is not coming
                # from an alternative settings. This allows support for
                # custom config files through populator without adding
                # a bunch of empty variable declarations in the run
                # script when using the standard config file
                self.input_filenames[inp_file_k] = ""

        sdos_values = { 'exe_path': self.l2_binary_filename,
                        'data_path': os.environ.get("L2_INPUT_PATH", None),
                        'src_path':  os.environ.get("L2_FP_SRC_PATH", None),
                        'algorithm_descriptor': 'L2 Full Physics',
                        'algorithm_maturity': 'operational',
                        'comments': 'Generated using L2 Populator' }

        self.logger.info("Writing input config values to: %s" %
                         self.input_config_filename)
        self.set_input_config_values(config_filename,
                                     self.config_sounding_id_section,
                                     self.input_config_filename,
                                     list(self.input_filenames.values()), **sdos_values)

        if self.l2_binary_filename:
            self.logger.info("Creating run scripts")
            self.create_run_scripts(config_filename,
                                    self.config_sounding_id_section)

        return True
       
 
