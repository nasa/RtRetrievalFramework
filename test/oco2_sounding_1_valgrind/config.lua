require "oco_base_config"

config = OcoBaseConfig:new()

-- For now, suppress use of EOFs. These are in bit of flux, and there is
-- no reason to change the output from this test each time we change the
-- EOFs.
config.fm.instrument.instrument_correction.ic_nadir = {}
config.fm.instrument.instrument_correction.ic_glint = {}
config.fm.instrument.instrument_correction.ic_target = {}

--- Newer OCO-2 configuration uses new fields in met file. But the old
--- test data doesn't have these fields. So use the older merra climatology
--- for these tests
config.fm.atmosphere.aerosol.creator = ConfigCommon.merra_aerosol_creator

config:do_config()
