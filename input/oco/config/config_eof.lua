
require "oco_base_config"

base_dir = ConfigCommon.local_dir()
config = OcoBaseConfig:new()

--- We look for the file in the same directory as this config lua file.
--- This is just a convenience, if it makes more sense you can just put
--- a full path in here and place the static input file wherever you like.
--- Note that this file has both the spectral window update and the EOF data
config.static_file = base_dir .. "/l2_oco_static_input_eof.h5"
config.fm.instrument.instrument_correction.ic = {"eof_1", "eof_2"}
config.fm.instrument.instrument_correction.eof_1 = {
   hdf_group = "Instrument/EmpiricalOrthogonalFunction",
   apriori = ConfigCommon.hdf_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
   covariance = ConfigCommon.hdf_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
   order = 1,
   by_pixel = true,
   creator = ConfigCommon.empirical_orthogonal_function,
   retrieve_bands = { true, true, true },
}
config.fm.instrument.instrument_correction.eof_2 = {
   hdf_group = "Instrument/EmpiricalOrthogonalFunction",
   apriori = ConfigCommon.hdf_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 2 - 1),
   covariance = ConfigCommon.hdf_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 2 - 1),
   order = 2,
   by_pixel = true,
   creator = ConfigCommon.empirical_orthogonal_function,
   retrieve_bands = { true, true, true },
}

config:do_config()
