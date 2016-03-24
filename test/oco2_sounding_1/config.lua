require "oco_base_config"

config = OcoBaseConfig:new()

-- For now, suppress use of EOFs. These are in bit of flux, and there is
-- no reason to change the output from this test each time we change the
-- EOFs.
config.fm.instrument.instrument_correction.ic_nadir = {}
config.fm.instrument.instrument_correction.ic_glint = {}
config.fm.instrument.instrument_correction.ic_target = {}

config:do_config()
