------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, but have Rayleigh only with
--- no Aerosol.
---
--- Also use LRad, but not LSI
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new()

--- Use LRAD only
config.fm.rt.nstream = 4
config.fm.rt.creator = ConfigCommon.radiative_transfer_lrad

--- Use Rayleigh only.
config.fm.atmosphere.aerosol.creator = ConfigCommon.rayleigh_only

config:do_config()
