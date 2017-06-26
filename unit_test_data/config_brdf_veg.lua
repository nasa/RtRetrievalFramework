------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, but with BreonVeg ground
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new() 

config.fm.atmosphere.ground.creator = ConfigCommon.ground_brdf_veg

config:do_config()
