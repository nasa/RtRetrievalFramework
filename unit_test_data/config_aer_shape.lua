------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, but with Coxmunk ground
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new()

aer_cfg = config.fm.atmosphere.aerosol
for idx, name in ipairs(aer_cfg.aerosols) do
   aer_cfg[name].apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log")
   aer_cfg[name].covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log")
   aer_cfg[name].creator = ConfigCommon.aerosol_log_shape_gaussian
end

config:do_config()
