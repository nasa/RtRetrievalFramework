------------------------------------------------------------
--- Sets up deviations from the default behavior for 
--- running orbit simulator r60 data
------------------------------------------------------------

require "oco_base_config"

config = OcoBaseConfig:new()

-- Use ECMWF files provided by CSU
--config.fm.input.met.creator = OcoConfig.oco_meteorology

-- ABSCO paths and filenames to match those used by the sim
config.fm.atmosphere.absorber.CO2.absco = "v4.0.2/lowres/co2_v4.0.2-lowres.hdf"
config.fm.atmosphere.absorber.H2O.absco = "v4.0.2/lowres/h2o_v4.0.2-lowres.hdf"
config.fm.atmosphere.absorber.O2.absco = "v4.0.2/lowres/o2_v4.0.2-lowres.hdf"

-- Match the same ILS extents as used in the simulator
config.fm.instrument.ils_half_width = { 
    DoubleWithUnit(1.50e-04, "um"), DoubleWithUnit(5.00e-04, "um"), DoubleWithUnit(5.00e-04, "um"), 
}

-- Old CO2 scale factor
old_aprior_func = config.fm.atmosphere.absorber.CO2.apriori
function config:scaled_apriori()
   local v = old_aprior_func(self)
   for i=0,v:rows() - 1 do
     v:set(i, v(i) * 0.9873)
   end
   return v
end
config.fm.atmosphere.absorber.CO2.apriori = config.scaled_apriori

-- Old O2 value
function config:orbit_sim_o2()
   local const_arr = Blitz_double_array_1d(self.config.pressure:max_number_level())
   const_arr:set(Range(), 0.2095)
   return const_arr
end
config.fm.atmosphere.absorber.O2.apriori = config.orbit_sim_o2

-- Use old solar model
static_solar_file = config_common_dir .. "/../input/old_orbit_sim_solar_model.h5"

config:do_config()
