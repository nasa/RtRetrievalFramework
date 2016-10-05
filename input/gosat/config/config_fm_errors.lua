require "gosat_base_config"

config = GosatBaseConfig:new()

config.absco_path = "/groups/algorithm/l2_fp/absco/"
config.fm.atmosphere.absorber.CO2.absco = "spectroscopy_errors/co2_v20120311.hdf"
config.fm.atmosphere.absorber.H2O.absco = "spectroscopy_errors/h2o_v20120311.hdf"
config.fm.atmosphere.absorber.O2.absco = "spectroscopy_errors/o2_v20120921.hdf"

config:do_config()
