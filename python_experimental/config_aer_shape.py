from full_physics import *
import numpy as np
import aerosol_shape

import os
import sys

###########

# In order to pass these to Lua, we must define the values here
# so they will exist in memory for the duration of the program
# execution.
GAUSSIAN_COEFFS = { 'Ice':     np.array([-3.28341, 0.3, 0.04]),
                    'Kahn_2b': np.array([-3.28341, 1.0, 0.2 ]),
                    'Kahn_3b': np.array([-3.28341, 1.0, 0.2]),
                    'Water':   np.array([-3.28341, 0.75, 0.1]),
                    }

GAUSS_COEFF_VAR = { 'Ice':     np.diag(np.square([1.5, 0.2, 0.01])),
                    'Kahn_2b': np.diag(np.square([1.5, 0.4, 0.01])),
                    'Kahn_3b': np.diag(np.square([1.5, 0.4, 0.01])),
                    'Water':   np.diag(np.square([1.5, 0.4, 0.01])),
                    }

###########

# Use the dynamic_config.lua as in a production run if we can find it, otherwise
# use unit test data for testing
prod_config = "./config/dynamic_config.lua"
if not os.path.exists(prod_config):
    small_set_dir = os.path.join(os.path.dirname(sys.argv[0]), "../unit_test_data/tccon_small_set")
    base_config = os.path.join(small_set_dir, "config/small_set_base.lua")
else:
    base_config = prod_config

ls = LuaState.load_file(os.path.realpath(base_config))
lg = ls.globals()
lg.config = lg.DynamicConfig.new(lg.DynamicConfig)
lua_config = ls.globals().config

if not os.path.exists(prod_config):
    # Use a specific sounding id from the small set, if not defined,
    # convenience during testing
    lua_config.sid_string = "20090827005603"

# Set up Lua configuration to use Gaussian shape model and data
def aerosol_gauss_apriori(c,aer_name):
    return GAUSSIAN_COEFFS[aer_name]

def aerosol_gauss_covariance(c,aer_name):
    return GAUSS_COEFF_VAR[aer_name]

class CreatorAerosolShape:
    def __init__(self, lg):
        self.c = lg.CreatorAerosol.new(lg.CreatorAerosol)
        self.c.extinction = self.extinction

    def extinction(self, c):
        return aerosol_shape.AerosolShapeGaussian(c.config.pressure, c.retrieval_flag(c), c.apriori_v(c), c.name)

aerosol_cfg = lua_config.fm.atmosphere.aerosol
for aer_name in aerosol_cfg.aerosols:
    aerosol_cfg[aer_name].apriori = aerosol_gauss_apriori
    aerosol_cfg[aer_name].covariance = aerosol_gauss_covariance
    aerosol_cfg[aer_name].creator = CreatorAerosolShape(lg).c
    
# Instruct Lua to create everything
lua_config.do_config(lua_config)
