from __future__ import absolute_import
from .acos_populator import AcosPopulator, PopulatorBase
import os

class GosatPopulator(AcosPopulator):
    '''This is a populator for Gosat'''
    def __init__(self, **user_settings):
        # Default value for this. Populator may have this overriden
        self.l2_config_filename = \
            os.path.join(os.environ.get('L2_INPUT_PATH', ''),
                         "gosat/config/config.lua")
        config_sounding_id_section = "input/GOSATFullPhysics/SoundingIds"

        config_input_keywords = \
        { 'spectrum_file': "input/InputProductFiles/L1BFile",
          'met_file':      "input/InputProductFiles/ResampledMetFile",
          'rrv_file':      "input/InputProductFiles/RRVFile",
        }
        AcosPopulator.__init__(self, config_sounding_id_section,
                               config_input_keywords, **user_settings)

    @property
    def config_file_obj_re(self):
        return r'^GOSATFullPhysics'

    @property
    def config_template_file(self):
       return os.path.dirname(__file__) + "/template/gosat_config.tmpl"


# Register class with PopulatorBase so it is known to populate
        
PopulatorBase.populator_list["gosat"] = GosatPopulator
        
        
