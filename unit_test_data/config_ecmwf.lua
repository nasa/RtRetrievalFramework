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
config.met_file = "in/ecmwf.h5"

--- Use HDF Level 1 file rather than ASCII
config.fm.input = {
    creator = ConfigCommon.l1b_met_input,
    l1b = {
        creator = AcosConfig.level1b_hdf,
        noise = {
            creator = ConfigCommon.noise_ascii_array
        },
    },
    met = {
        creator = AcosConfig.acos_met,
    },
}

config.fm.atmosphere.pressure.apriori = ConfigCommon.met_pressure
config.fm.atmosphere.temperature.levels.apriori = ConfigCommon.met_temperature_fixed_pressure
config.fm.atmosphere.absorber.H2O.apriori = ConfigCommon.met_h2o_vmr_fixed_pressure

config:do_config()
