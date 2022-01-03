------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, but with Scaled Coxmunk ground
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new() 

config.fm.atmosphere.ground.creator = ConfigCommon.ground_coxmunk_scaled

config:do_config()
