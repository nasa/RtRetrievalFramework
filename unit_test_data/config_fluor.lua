------------------------------------------------------------
--- Adds the use of FluorescenceEffect to the retrieval
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new()

function fluorescence_ap()
   local f_coeffs = Blitz_double_array_1d(2)
   f_coeffs:set(0, -1.35039e-09)
   f_coeffs:set(1, 0.0016)
   return f_coeffs
end

function fluorescence_cov()
   local f_cov = Blitz_double_array_2d(2,2)
   f_cov:set(Range.all(), Range.all(), 0)
   f_cov:set(0, 0, 0.02 * 0.02)
   f_cov:set(1, 1, 7e-4 * 7e-4)
   return f_cov
end

function reference_point_func()
   return DoubleWithUnit(0.755, "micron")
end

config.fm.spectrum_effect.fluorescence = {
   apriori = fluorescence_ap,
   covariance = fluorescence_cov,
   creator = ConfigCommon.fluorescence_effect,
   reference_point = reference_point_func,
   retrieved = true,
}
table.insert(config.fm.spectrum_effect.speceff, "fluorescence")

config:do_config()
