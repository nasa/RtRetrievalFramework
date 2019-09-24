require "oco_base_config"

config = OcoBaseConfig:new()

config.fm.spec_win.bad_sample_mask = nil

-- For now, suppress use of EOFs. These are in bit of flux, and there is
-- no reason to change the output from this test each time we change the
-- EOFs.
config.fm.instrument.instrument_correction.ic_nadir = {"zero_offset_waveform"}
config.fm.instrument.instrument_correction.ic_glint = {"zero_offset_waveform"}
config.fm.instrument.instrument_correction.ic_target = {"zero_offset_waveform"}

--- Newer OCO-2 configuration uses new fields in met file. But the old
--- test data doesn't have these fields. So use the older merra climatology
--- for these tests
config.fm.atmosphere.aerosol.creator = ConfigCommon.merra_aerosol_creator

--- Newer OCO-2 configuration gets the CO2 prior from a L2CPr file. But
--- old test data doesn't have this file, use the older
---- reference_co2_apriori_met_apriori
config.fm.atmosphere.absorber.CO2.apriori = ConfigCommon.reference_co2_apriori_met_apriori

config:do_config()
