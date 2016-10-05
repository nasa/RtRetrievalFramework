from full_physics import *
import aerosol_log_10

# Get Defaults
ls = LuaState.load_file("../unit_test_data/full_run/base_config.lua")

lg = ls.globals()
lg.config = lg.BaseConfig.new(lg.BaseConfig)
lua_config = ls.globals().config

# Local modifications
lua_config.create_aerosol = aerosol_log_10.create_aerosol_log_10

# Now create everything
lua_config.do_config(lua_config)


