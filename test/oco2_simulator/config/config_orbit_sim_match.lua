------------------------------------------------------------
--- This matches a specific orbit simulator run. For this
--- run we have Rayleigh scattering only. We base all our
--- initial guess to match what was actually used the OCO
--- simulator run, so if we exactly match the OCO simulator
--- in Level 2 this will give a 0 residual with the initial
--- state vector.
------------------------------------------------------------

require "oco_baseline_config"

config = OcoBaselineConfig:new()

-- Turn off aerosol retrieval and scattering calculations
config.fm.atmosphere.aerosol.creator = ConfigCommon.rayleigh_only

config.fm.instrument.dispersion.retrieved=true

------------------------------------------------------------
--- Use meteorology file for ECMWF
------------------------------------------------------------

config.fm.input.met.creator = OcoConfig.oco_meteorology

------------------------------------------------------------
--- Use apriori that matches what was used in the scene file
------------------------------------------------------------

config.fm.atmosphere.absorber.CO2.apriori = OcoConfig.co2_apriori_from_scene


-------------------------------------------------------------
---- Orbit simulator used a much more narrow ILS than normal OCO.
---- This was later changed in the simulator, but to
---- match the data we currently have we need to match what
---- was done.
-------------------------------------------------------------
config.fm.instrument.ils_half_width = { DoubleWithUnit(1.5e-4, "um"), 
                                        DoubleWithUnit(5.0e-4, "um"),
                                        DoubleWithUnit(5.0e-4, "um") }

------------------------------------------------------------
--- Set lambertian albedo to the value we know was set in the
--- simulator.
------------------------------------------------------------

config.fm.atmosphere.ground.lambertian.apriori = 
   ConfigCommon.fixed_albedo({0.3, 0.3, 0.3}, {0.0, 0.0, 0.0})

-- No EOFs
config.fm.instrument.instrument_correction.ic_nadir = {}
config.fm.instrument.instrument_correction.ic_glint = {}
config.fm.instrument.instrument_correction.ic_target = {}

config:do_config()
