from __future__ import absolute_import
from .acos_populator import AcosPopulator, PopulatorBase
import os

class OcoUplookingPopulator(AcosPopulator):
    '''This is a populator for Oco uplooking'''
    def __init__(self, **user_settings):
        # Default value for this. Populator may have this overriden
        self.l2_config_filename = \
            os.path.join(os.environ.get('L2_INPUT_PATH', ''),
                         "oco/config/config_uplooking.lua")
        config_sounding_id_section = "input/OCOFullPhysics/SoundingIds"

        config_input_keywords = \
        { 'spectrum_file': "input/InputProductFiles/L1BFile",
          'atmosphere_file':    "input/InputProductFiles/AtmosphereFile" }
        AcosPopulator.__init__(self, config_sounding_id_section,
                               config_input_keywords, **user_settings)

    @property
    def config_file_obj_re(self):
        return r'^OCOUplookingFullPhysics'

    @property
    def config_template_file(self):
       return os.path.dirname(__file__) + "/template/oco_uplooking_config.tmpl"
        

# Register class with PopulatorBase so it is known to populate
        
PopulatorBase.populator_list["oco_uplooking"] = OcoUplookingPopulator
