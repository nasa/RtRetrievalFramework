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

-- Normal static input file, except we have additional aerosol types in it
config = GosatBaseConfig:new {
   static_file = gosat_base_config_dir .. "/../input/l2_gosat_static_input_odell.h5",
}
-- Chris modification

-- Different aerosol list
config.fm.atmosphere.aerosol.aerosols = {"ic_060", "du_1_3", "su_13_1"}
config.fm.atmosphere.aerosol.ic_060 = {
   creator = ConfigCommon.aerosol_log_shape_gaussian,
   apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
   covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
   property = ConfigCommon.hdf_aerosol_property("Aerosol"),
}
config.fm.atmosphere.aerosol.du_1_3 = {
   creator = ConfigCommon.aerosol_log_shape_gaussian,
   apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
   covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
   property = ConfigCommon.hdf_aerosol_property("Aerosol"),
}
config.fm.atmosphere.aerosol.su_13_1 = {
   creator = ConfigCommon.aerosol_log_shape_gaussian,
   apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
   covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
   property = ConfigCommon.hdf_aerosol_property("Aerosol"),
}

-- Turn on Jacobians in output
config.write_jacobian = true

config:do_config()
