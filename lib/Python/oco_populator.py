from __future__ import absolute_import
from .acos_populator import AcosPopulator, PopulatorBase
import os

class OcoPopulator(AcosPopulator):
    '''This is a populator for Oco'''
    def __init__(self, **user_settings):
        # Default value for this. Populator may have this overriden
        self.l2_config_filename = \
            os.path.join(os.environ.get('L2_INPUT_PATH', ''),
                         "oco/config/config.lua")
        config_sounding_id_section = "input/OCOFullPhysics/SoundingIds"

        config_input_keywords = \
        { 'spectrum_file':      "input/InputProductFiles/L1BFile",
          'met_file':           "input/InputProductFiles/ResampledMetFile",
          'co2_pr_file':        "input/InputProductFiles/CO2PriorFile",
          'imap_file':          "input/InputProductFiles/IMAPFile",
          'aband_file':         "input/InputProductFiles/ABandFile",
          'rrv_file':           "input/InputProductFiles/RRVFile",
          'input_file_mapping': "input/InputProductFiles/InputFileMapping",
          }
        AcosPopulator.__init__(self, config_sounding_id_section,
                               config_input_keywords,
                               **user_settings)

    @property
    def config_file_obj_re(self):
        return r'^OCOFullPhysics'

    @property
    def config_template_file(self):
       return os.path.dirname(__file__) + "/template/oco_config.tmpl"

class Oco3Populator(AcosPopulator):
    '''This is a populator for Oco3'''
    def __init__(self, **user_settings):
        # Default value for this. Populator may have this overriden
        self.l2_config_filename = \
            os.path.join(os.environ.get('L2_INPUT_PATH', ''),
                         "oco/config/config_oco3.lua")
        config_sounding_id_section = "input/OCO3FullPhysics/SoundingIds"

        config_input_keywords = \
        { 'spectrum_file':      "input/InputProductFiles/L1BFile",
          'met_file':           "input/InputProductFiles/ResampledMetFile",
          'co2_pr_file':        "input/InputProductFiles/CO2PriorFile",
          'imap_file':          "input/InputProductFiles/IMAPFile",
          'aband_file':         "input/InputProductFiles/ABandFile",
          'rrv_file':           "input/InputProductFiles/RRVFile",
          'input_file_mapping': "input/InputProductFiles/InputFileMapping",
          }
        AcosPopulator.__init__(self, config_sounding_id_section,
                               config_input_keywords,
                               **user_settings)

    @property
    def config_file_obj_re(self):
        return r'^OCO3FullPhysics'

    @property
    def config_template_file(self):
       return os.path.dirname(__file__) + "/template/oco3_config.tmpl"
   
# Register class with PopulatorBase so it is known to populate
        
PopulatorBase.populator_list["oco"] = OcoPopulator
PopulatorBase.populator_list["oco3"] = Oco3Populator
        
