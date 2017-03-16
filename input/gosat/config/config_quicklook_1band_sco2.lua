------------------------------------------------------------
--- This is the configuration used in the "quick look" test
--- (the 25,000 soundings test). 
---
--- If these changes get approved, they will need to be moved
--- into the production configuration file (e.g., gosat_base_config),
--- but we keep these separate while testing
------------------------------------------------------------

--- Test B3.4.00_quick_look_testing_9, turn off fluorescence for M gain

--- Test B3.4.00_quick_look_testing_7, try 2 different ILS updates for band
---  A (ILS_1)
--- Test B3.4.00_quick_look_testing_6, use 2nd order polynomial for albedo
--- (change turned back off)
--- Test B3.4.00_quick_look_testing_5, have input file that changes pressure
--- constraint, turn off EOF
--- Test B3.4.00_quick_look_testing_4, change version of ABSCO tables with
--- new pressure gridding
--- Test B3.4.00_quick_look_testing_3, have input file that changes WCO2
--- spectral range to match TCCON

require "gosat_base_config"
require "single_band_support"

config = GosatBaseConfig:new {
   static_file = gosat_base_config_dir .. "/../input/l2_gosat_static_input_quicklook.h5",
}
--- Test 4b, new ABSCO tables
config.fm.atmosphere.absorber.CO2.absco = "v4.2.0_beta_narrow/co2_v4.2.0_with_ctm_sco2_rescaled_by_0.99.hdf"
config.fm.atmosphere.absorber.H2O.absco = "v4.2.0_beta_narrow/h2o_v4.2.0.hdf"
config.fm.atmosphere.absorber.O2.absco = "v4.2.0_beta_narrow/o2_v4.2.0_drouin_rescaled_by_1.0125.hdf"
-- Test 5, Turn off EOF
config.fm.instrument.instrument_correction.ic = { "zero_offset_waveform" }
--- Test 6, 2nd order polynomial for albedo
--- (Turned back off)
-- config.fm.atmosphere.ground.lambertian.apriori =
--   AcosConfig.gosat_albedo_from_radiance(2)
--- Test 7a theoretical ILS_1
--- Settled on theoretical ILS
config.fm.instrument.ils_func.hdf_group_name = "Instrument/ILS_theoretical"
--- Test 7b latest JAXA update of ILS_1
-- config.fm.instrument.ils_func.hdf_group_name = "Instrument/ILS_jaxa_update"

-- Test 9, turn off fluorescence for M gain
-- config.fm.spectrum_effect.creator = AcosConfig.spectrum_effect_list_h_and_m
-- config.fm.spectrum_effect.speceff_h_gain = { "solar_model", "fluorescence" }
-- config.fm.spectrum_effect.speceff_m_gain = { "solar_model",  }

-- Set up retrieval as suggested by Christian
--- Use a scaled version of the CO2 retrieval

config.fm.atmosphere.absorber.CO2.vmr_profile = ConfigCommon.tccon_co2_apriori_met
config.fm.atmosphere.absorber.CO2.scale_apriori = 1.0
config.fm.atmosphere.absorber.CO2.scale_cov = 0.1
config.fm.atmosphere.absorber.CO2.creator = ConfigCommon.vmr_level_scaled

-- Rayleigh only
config.fm.atmosphere.aerosol.creator = ConfigCommon.rayleigh_only

-- Turn off zero level offset and fluorescence
config.fm.instrument.instrument_correction.ic = {  }
config.fm.spectrum_effect.speceff = { "solar_model" }

-- Hold pressure fixed
config.fm.atmosphere.pressure.retrieved = false

-- Use a difference lambertian term over the ocean (we use the default values
-- over land)
config.fm.atmosphere.ground.coxmunk_lambertian.apriori = ConfigCommon.hdf_apriori_i("Ground/Coxmunk_Albedo_single_band")
config.fm.atmosphere.ground.coxmunk_lambertian.covariance = ConfigCommon.hdf_covariance_i("Ground/Coxmunk_Albedo_single_band")

-- Do a single band retrieval
config.which_spectrometers = "SCO2"
init_single_band_support(config)

config:do_config()
