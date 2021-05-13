require "oco3_base_config"

config = Oco3BaseConfig:new()

config.diagnostic = true

config.fm.atmosphere.aerosol.creator = ConfigCommon.merra_aerosol_creator
config.fm.atmosphere.absorber.CO2.apriori = ConfigCommon.reference_co2_apriori_met_apriori

config:do_config()
