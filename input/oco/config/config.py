import os
from full_physics import *

# Get Defaults
config_dir = os.path.dirname(os.path.realpath(__file__))
ls = LuaState.load_file(str(os.path.join(config_dir, "base_config.lua")))

lg = ls.globals()
lg.config = lg.BaseConfig.new(lg.BaseConfig)
lua_config = ls.globals().config

# Now create everything
lua_config.do_config(lua_config)
