------------------------------------------------------------
-- Makes the  necessary changes to the nominal OCO
-- configuration to support uplooking data retrievals.
------------------------------------------------------------

dofile "ascii_atm_base_config.lua"

config = AsciiAtmBaseConfig:new()

config.do_retrieval = true
config.sid_string = "2007032016155813"
config.spectrum_file = "../input/l1b_simulated.h5"
config.atmosphere_file = "../input/atmosphere_00999_20.dat"
config.albedo_file = "../input/albedo_00999.dat"

-- Extend the pressure levels bottom level to give more room for surface
-- pressure movement
orig_press_levs = config.fm.atmosphere.pressure.pressure_levels
function extended_pressure_level(self)
   plev = orig_press_levs(self)
   plev:set(plev:rows()-1, 105000)
   return plev
end
config.fm.atmosphere.pressure.pressure_levels = extended_pressure_level

config:do_config()
