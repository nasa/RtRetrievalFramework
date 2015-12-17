------------------------------------------------------------
-- Adds necessary changes for using radiance scaling 
-- instead of coxmunk + lamberitan
------------------------------------------------------------

--- *********************************************************
---  Note, this is now the default configuration, and much of
---  This was moved to acs_config.lua. We'll leave this in 
---  place for now, but we shouldn't need to use this anymore
--- *********************************************************

require "gosat_base_config"

config = GosatBaseConfig:new()

table.insert(config.fm.instrument.instrument_correction.ic, "radiance_scaling")

config:do_config()
