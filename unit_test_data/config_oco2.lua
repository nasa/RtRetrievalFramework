------------------------------------------------------------
--- Note that the use of base_config.lua is entirely optional.
--- There is a small number of global variables that you
--- need to create, you can do this however you want.
---
--- The variables are:
--- logger - A LogImp
--- forward_model - A ForwardModel
--- solver - A ConnorSolver
--- initial_guess - A InitialGuess
--- number_pressure_level - Integer giving the number of pressure 
---    levels. This is used to size the output file, so it should
---    be the maximim pressure
--- number_aerosol - Integer giving number of Aerosol particles. 
---    Can be 0. This is used to size the output file.
--- number_band - Integer giving number of spectral bands.
--- iteration_output - Boolean. True if we should write out iteration
---   output
--- register_output - A VectorRegisterOutput giving the
---   list of output that should be generated. This list can empty if
---   no output is desired.
------------------------------------------------------------

require "oco_base_config"

config = OcoBaseConfig:new()

config.sid_string = "2010090900133834"
config.cloud_file = ""
config.spectrum_file = "oco2_L1bScND_80008a_111017225030d_spliced.h5"
config.met_file = "oco2_ECMWFND_80008a_111018214952d_spliced.h5"
config.static_file = "../input/oco/input/l2_oco_static_input.h5"

config.fm.atmosphere.absorber.CO2.absco = "v3.3.0/lowres/co2_v3.3.0-lowres.hdf"
config.fm.atmosphere.absorber.H2O.absco = "v3.3.0/lowres/h2o_v3.3.0-lowres.hdf"
config.fm.atmosphere.absorber.O2.absco  = "v3.3.0/lowres/o2_v3.3.0-lowres.hdf"

--- Only use 3 EOFs for now
config.fm.instrument.instrument_correction.ic_nadir = {"eof_glint_1", "eof_glint_2","eof_glint_3"}
config.fm.instrument.instrument_correction.ic_glint = {"eof_glint_1", "eof_glint_2","eof_glint_3"}
config.fm.instrument.instrument_correction.ic_target = {"eof_glint_1", "eof_glint_2","eof_glint_3"}

--- Newer OCO-2 configuration uses new fields in met file. But the old
--- test data doesn't have these fields. So use the older merra climatology
--- for these tests
config.fm.atmosphere.aerosol.creator = ConfigCommon.merra_aerosol_creator

--- Newer OCO-2 configuration gets the CO2 prior from a L2CPr file. But
--- old test data doesn't have this file, use the older
---- reference_co2_apriori_met_apriori
config.fm.atmosphere.absorber.CO2.apriori = ConfigCommon.reference_co2_apriori_met_apriori

config:do_config()
