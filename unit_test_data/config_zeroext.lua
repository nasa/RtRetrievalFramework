------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, but read a zero extinction 
--- for each aerosol from an ascii file. Also use LRad, but
--- not LSI
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new()

--- Use zero extent.
for i,a in ipairs(config.fm.atmosphere.aerosol.aerosols) do
   config.fm.atmosphere.aerosol[a].apriori = ConfigCommon.ascii_aerosol_apriori("./in/aerosol_zeroext_20.dat")
   config.fm.atmosphere.aerosol[a].creator = ConfigCommon.aerosol_linear_profile
end

--- Use LRAD only
config.fm.rt.nstream = 4
config.fm.rt.creator = ConfigCommon.radiative_transfer_lrad

config:do_config()
