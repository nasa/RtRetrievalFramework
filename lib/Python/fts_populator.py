from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from .populator_base import PopulatorBase
from .l2_input import L2InputFile
import os
import numpy
from full_physics.oco_matrix import OcoMatrix

class FtsPopulator(PopulatorBase):
    '''This is a populator that is used for FTS.'''

    def __init__(self, **user_settings):
        # Default value for this. Populator may have this overriden
       self.l2_config_filename = \
            os.path.join(os.environ.get('L2_INPUT_PATH', ''),
                         "fts/config/config.lua")
       PopulatorBase.__init__(self, **user_settings)
       # We don't have a L1B1 file, so tell PopulatorBase not to try
       # to aggregate this.
       self.have_l1b1 = False
       self.pout_col_idx   = 27
       self.convert_factor = 1e2

    @property
    def spectrum_a_list_filename(self):
        return self.processing_dir + "/spectrum_a.list"

    @property
    def spectrum_b_list_filename(self):
        return self.processing_dir + "/spectrum_b.list"

    @property
    def variable_exports_list(self):
        return { "runlog_file": self.runlog_filename,
                 "atmosphere_file": self.atmosphere_filename,
                 }

    @property
    def additional_var_list(self):
        return { 'spectrum_a_file': self.spectrum_a_list_filename,
                 'spectrum_b_file': self.spectrum_b_list_filename,
                 }

    @property
    def config_file_obj_re(self):
        return r'^FTSFullPhysics'

    @property
    def config_template_file(self):
       return os.path.dirname(__file__) + "/template/fts_config.tmpl"

    # This doesn't appear to be used anymore, but keep the code here in
    # case we need it in the future.
    def create_psurf_apriori_file(self, spectrum_filename,
                                  psurf_output_filename):
        base_spec_name = os.path.basename(spectrum_filename)

        # Use grep because its faster than doing it outself
        grep_cmd = "grep -E " + base_spec_name + " " + self.runlog_filename
        matched_line = os.popen(grep_cmd).readline()

        if matched_line == None or len(matched_line) == 0:
            raise IOError('Could not find spectrum name: %s in run log file: %s' % (base_spec_name, self.runlog_filename))

        try:
            matched_columns = matched_line.split()
            psurf_val = float(matched_columns[self.pout_col_idx]) * self.convert_factor
        except:
            raise ValueError('Failed to parse psurf value from: "%s" from runlog line: %s' % (matched_columns[self.pout_col_idx], matched_line))
   
        out_obj = OcoMatrix()
    
        out_obj.data = numpy.zeros((1,1), dtype=float)
        out_obj.data[0,0] = psurf_val

        out_obj.file_id = 'psurf value extracted for spectrum named: %s from runlog file: %s' % (base_spec_name, self.runlog_filename)
        out_obj.labels = ['PSURF']
        
        out_obj.write(psurf_output_filename)

    def read_map_values_file(self, map_filename, section=None):
        '''Not sure what this does'''

        # Use old L2_Input syntax
        section = section.replace("/", "->")

        map_obj = L2InputFile(map_filename)

        if section != None:
            self.logger.debug('Reading map from section %s file: %s' % (section, map_filename))
            section_obj = map_obj.get_section(section)
        else:
            self.logger.debug('Reading id list from file: %s' % id_list_file)
            section_obj = map_obj.rootNode

        if section_obj == None or len(section_obj) == 0:
            raise IOError('Could not find section %s in file: %s' % (setion, map_filename))

        map_values = {}
        for currfilesect in section_obj:
            for sectkeyname in currfilesect.get_all_keyword_names():
                sectkeyval = currfilesect.get_keyword_value(sectkeyname)
                map_values[str(sectkeyname)] = str(sectkeyval)

        return map_values

    def populate(self, config_filename):
        '''This is the function that takes the configuration file and
        user settings, and use this to generate the run scripts.'''
        self.processing_dir = os.path.dirname(os.path.abspath(config_filename))
        config_obs_id_section = "input/FTSFullPhysics/FTSObsIds"
        config_runlog_file_section = "input/InputProductFiles/RunlogFile"
        config_atmosphere_file_section = "input/InputProductFiles/AtmosphereFile"
        config_spectrum_files_section = "input/InputProductFiles/SpectrumFiles"

        # Check if the configuration can be processed, else signal
        # that we can not proceed
        self.logger.info("Checking if processable")
        if not self.is_processable(config_filename, config_obs_id_section):
            return False

        # Create common necessary files and directories
        self.logger.info("Initializing processing dir")
        self.init_processing_dir()

        # Get observation ids
        obs_ids = self.read_id_list_file(config_filename, config_obs_id_section)

        # Create lists of spectrum files for run script
        self.logger.info("Creating spectrum filename lists")
        spectrum_files_map = self.read_map_values_file(config_filename, config_spectrum_files_section)

        spec_a_out = open(self.spectrum_a_list_filename, "w")
        spec_b_out = open(self.spectrum_b_list_filename, "w")
        for curr_id in obs_ids:
            print(spectrum_files_map[curr_id + "1"], file=spec_a_out)
            print(spectrum_files_map[curr_id + "2"], file=spec_b_out)
        spec_a_out.close()
        spec_b_out.close()

        # File paths needed by configuration
        self.logger.info("Creating surface pressure apriori files")

        self.runlog_filename = self.get_config_keyword_value(config_filename, config_runlog_file_section)
        self.atmosphere_filename = self.get_config_keyword_value(config_filename, config_atmosphere_file_section)

        if self.l2_binary_filename:
            self.logger.info("Creating run scripts")
            self.create_run_scripts(config_filename, config_obs_id_section)

        return True

# Register class with PopulatorBase so it is known to populate
        
PopulatorBase.populator_list["fts"] = FtsPopulator
