------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, but with Coxmunk ground
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new()

config.fm.atmosphere.ground.creator = ConfigCommon.ground_coxmunk

config:do_config()
