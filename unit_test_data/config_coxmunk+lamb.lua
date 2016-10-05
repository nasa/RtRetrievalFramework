------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, but with Coxmunk+Lamb ground
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new() 

config.fm.atmosphere.ground.creator = ConfigCommon.ground_coxmunk_plus_lamb

config:do_config()
