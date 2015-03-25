------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, but with BreonVeg ground
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new() 

config.fm.atmosphere.ground.creator = ConfigCommon.ground_breon_veg

config:do_config()
