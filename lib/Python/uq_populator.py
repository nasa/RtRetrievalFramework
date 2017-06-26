from __future__ import absolute_import
from .acos_populator import AcosPopulator, PopulatorBase
import os

class UqPopulator(AcosPopulator):
    '''This is a populator for Uncertainty Quantification testing'''
    def __init__(self, **user_settings):
        # Default value for this. Populator may have this overriden
        self.l2_config_filename = \
            os.path.join(os.environ.get('L2_INPUT_PATH', ''),
                         "oco/config/config_uq.lua")
        config_sounding_id_section = "input/UqFullPhysics/SoundingIds"

        config_input_keywords = \
        { 'spectrum_file': "input/InputProductFiles/UqFile",
        }

        AcosPopulator.__init__(self, config_sounding_id_section,
                               config_input_keywords,
                               **user_settings)

        # Don't create an aggregation script that includes the L1B file
        self.have_l1b1 = False

    @property
    def config_file_obj_re(self):
        return r'^UqFullPhysics'

    @property
    def config_template_file(self):
       return os.path.dirname(__file__) + "/template/uq_config.tmpl"
        
# Register class with PopulatorBase so it is known to populate
        
PopulatorBase.populator_list["uq"] = UqPopulator
