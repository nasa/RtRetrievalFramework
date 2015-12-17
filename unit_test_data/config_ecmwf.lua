------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, with changes listed below
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new {

---- Input for HDF file
   noise_file = "in/l1b/spec/spectra.dat",
   spectrum_file = "in/sounding_id.h5",
   sid_string = "20091009203401",

--- No cloud to go with this.
   cloud_file = "",
}
--- Pressure from ECMWF
config.ecmwf_file = "in/ecmwf.h5"

--- Use HDF Level 1 file rather than ASCII
config.fm.l1b.creator = AcosConfig.level1b_hdf
config.fm.l1b.noise.creator = ConfigCommon.noise_ascii_array

config.fm.atmosphere.pressure.apriori = ConfigCommon.ecmwf_pressure
config.fm.atmosphere.temperature.levels.apriori = ConfigCommon.ecmwf_temperature
config.fm.atmosphere.absorber.H2O.apriori = ConfigCommon.ecmwf_h2o_vmr

config:do_config()
