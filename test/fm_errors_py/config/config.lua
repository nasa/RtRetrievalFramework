------------------------------------------------------------
-- Makes the  necessary changes to the nominal OCO2
-- configuration using a compact test set from the
-- 20120516 repackaged OCO2 orbit simulator data sample
------------------------------------------------------------

require "oco_base_config"

config = OcoBaseConfig:new()
config.sid_string = "2010090911454303"
config.spectrum_file = "../input/oco2_L1bScND_89012a_100909_Bxxxx_120315145001d_spliced.h5"
config.met_file = "../input/oco2_ECMWFND_89012a_100909_Bxxxx_120510145007d_spliced.h5"

-- Ensure that low stream = 1 (ie 2 full stream) for this test
config.fm.rt.lsi_constant.low_stream = 1 

-- Use dedicated 2stream if not already enabled
config.fm.rt.lsi_constant.dedicated_twostream = true

-- Use old ILS extents
config.fm.instrument.ils_half_width = { 
    DoubleWithUnit(1.2e-3, "um"), DoubleWithUnit(5.1e-3, "um"), DoubleWithUnit(8.4e-3, "um"), 
}

-- Use old PSURF covariance
function old_oco_psurf_cov_value(self) 
    local cov = Blitz_double_array_2d(1,1)
    cov:set(0, 0, 625)
    return cov
end
config.fm.atmosphere.pressure.covariance = old_oco_psurf_cov_value

--- No EOFs
config.fm.instrument.instrument_correction.ic_nadir = {}
config.fm.instrument.instrument_correction.ic_glint = {}
config.fm.instrument.instrument_correction.ic_target = {}

config:do_config()
