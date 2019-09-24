------------------------------------------------------------
-- Makes the  necessary changes to the nominal OCO2
-- configuration using a compact test set from the
-- 20120516 repackaged OCO2 orbit simulator data sample
------------------------------------------------------------

require "oco_base_config"

config = OcoBaseConfig:new()

config.do_retrieval = false
config.sid_string = "2010090911454303"
config.spectrum_file = "../input/oco2_L1bScND_89012a_100909_Bxxxx_120315145001d_spliced.h5"
config.met_file = "../input/oco2_ECMWFND_89012a_100909_Bxxxx_120510145007d_spliced.h5"

-- Ensure that low stream = 1 (ie 2 full stream) for this test
config.fm.rt.lsi_constant.low_stream = 1 

config.fm.atmosphere.absorber.CO2.absco = "v3.3.0/lowres/co2_v3.3.0-lowres.hdf"
config.fm.atmosphere.absorber.H2O.absco = "v3.3.0/lowres/h2o_v3.3.0-lowres.hdf"
config.fm.atmosphere.absorber.O2.absco = "v3.3.0/lowres/o2_v3.3.0-lowres.hdf"

-- Use old ILS extents
config.fm.instrument.ils_half_width = { 
    DoubleWithUnit(1.2e-3, "um"), DoubleWithUnit(5.1e-3, "um"), DoubleWithUnit(8.4e-3, "um"), 
}

-- No EOFs
config.fm.instrument.instrument_correction.ic_nadir = {}
config.fm.instrument.instrument_correction.ic_glint = {}
config.fm.instrument.instrument_correction.ic_target = {}

--- Newer OCO-2 configuration uses new fields in met file. But the old
--- test data doesn't have these fields. So use the older merra climatology
--- for these tests
config.fm.atmosphere.aerosol.creator = ConfigCommon.merra_aerosol_creator

--- Newer OCO-2 configuration gets the CO2 prior from a L2CPr file. But
--- old test data doesn't have this file, use the older
---- reference_co2_apriori_met_apriori
config.fm.atmosphere.absorber.CO2.apriori = ConfigCommon.reference_co2_apriori_met_apriori
