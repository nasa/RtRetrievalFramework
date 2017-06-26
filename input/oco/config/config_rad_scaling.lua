------------------------------------------------------------
-- Adds necessary changes for using radiance scaling 
-- instead of coxmunk + lamberitan
------------------------------------------------------------

require "oco_base_config"

config = OcoBaseConfig:new()

able.insert(config.fm.instrument.instrument_correction.ic, "radiance_scaling")

config:do_config()
