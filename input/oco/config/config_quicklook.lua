------------------------------------------------------------
--- This is the changes needed to run the oco2_simulator
--- tests
------------------------------------------------------------

require "oco_base_config"

config = OcoBaseConfig:new {
--- static_file = oco_base_config_dir .. "/../input/l2_oco_static_input_quicklook.h5",
}

-- == Configuration deduced from /scratch2/algorithm/OCO2_orbit_simulations/r69/02c/nadir_glint/20120618_170/OCO2_sim_NDa_20120618_170_r69CSUSim02c_020.log ==

-- ABSCO paths and filenames
config.fm.atmosphere.absorber.CO2.absco = "v4.1.1/lowres/co2_v4.1.1-lowres.hdf"
config.fm.atmosphere.absorber.H2O.absco = "v4.1.1/lowres/h2o_v4.1.1-lowres.hdf"
config.fm.atmosphere.absorber.O2.absco = "v4.1.1/lowres/o2_v4.1.1-lowres.hdf"
config.fm.atmosphere.absorber.CO2.table_scale = 1.0
config.fm.atmosphere.absorber.O2.table_scale = 1.0

-- Match the same ILS extents as used in the simulator
config.fm.instrument.ils_half_width = { 
    DoubleWithUnit(4.09e-04, "um"), DoubleWithUnit(1.08e-03, "um"), DoubleWithUnit(1.40e-03, "um"), 
}

-- == END ==

config:do_config()
