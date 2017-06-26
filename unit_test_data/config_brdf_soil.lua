------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, but with BreonSoil ground
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new() 

config.fm.atmosphere.ground.creator = ConfigCommon.ground_brdf_soil

config:do_config()
