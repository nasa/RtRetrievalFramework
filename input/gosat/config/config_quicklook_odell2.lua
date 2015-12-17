------------------------------------------------------------
--- This is the configuration used in the "quick look" test
--- (the 25,000 soundings test). 
---
--- If these changes get approved, they will need to be moved
--- into the production configuration file (e.g., gosat_base_config),
--- but we keep these separate while testing
------------------------------------------------------------

-- This is done off of B3.4

require "gosat_base_config"
require "single_band_support"

config = GosatBaseConfig:new()
-- Chris modification

-- Rayleigh only
config.fm.atmosphere.aerosol.creator = ConfigCommon.rayleigh_only

-- Exclude O2
config.which_spectrometers = "WCO2 SCO2"
init_single_band_support(config)

-- Hold pressure fixed
config.fm.atmosphere.pressure.retrieved = false

-- Turn off zero level offset and fluorescence
config.fm.instrument.instrument_correction.ic = {  }
config.fm.spectrum_effect.speceff_h_gain = { "solar_model"}

-- Use lambertian everywhere
config.fm.atmosphere.ground.creator = ConfigCommon.ground_lambertian

config:do_config()
