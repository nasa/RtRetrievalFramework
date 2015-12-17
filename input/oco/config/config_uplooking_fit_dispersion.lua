------------------------------------------------------------
-- Makes the  necessary changes to the nominal OCO
-- configuration to support uplooking data retrievals.
------------------------------------------------------------

dofile "base_uplooking_tvac2.lua"

config = tvac2_config:new()

config.diagnostic = true
config.solver_constant.max_iteration = 10
config.solver_constant.max_divergence = 2

-- One band at a time
spectral_window_one_band_only = Creator:new()
function spectral_window_one_band_only:create()
   orig_range = self.config:h():read_double_with_unit_3d("Spectral_Window/microwindow")
   -- "zero" out the other bands
   orig_range.data:set(0, 0, Range.all(), 0.0)
   --orig_range.data:set(1, 0, Range.all(), 0.0)
   orig_range.data:set(2, 0, Range.all(), 0.0)
   return SpectralWindowRange(orig_range)
end
--config.fm.spec_win.creator = spectral_window_one_band_only

-- Turn off atmospheric constituent retrievals
config.fm.atmosphere.pressure.retrieved = false
config.fm.atmosphere.temperature.retrieved = false
config.fm.atmosphere.absorber.H2O.retrieved = false

-- Make CO2 covariance so small it doesnt change
--config.fm.atmosphere.absorber.CO2.retrieved = false
config.fm.atmosphere.absorber.CO2.scale_cov = 1.0e-20

-- Start with a closer radiance scaling ap
function radiance_scaling_ap(self, i)
   ap = Blitz_double_array_1d(2)
   ap:set(Range.all(), 0)
   if i == 0 then
      ap:set(0, 0.04292775836)
   elseif i == 1 then
      ap:set(0, 0.05109568608)
   else
      ap:set(0, 0.05398459837)
   end
   return ap
end
--config.fm.instrument.instrument_correction.radiance_scaling.apriori = radiance_scaling_ap   
--config.fm.instrument.instrument_correction.radiance_scaling.retrieve_bands = { false, true, false }

-- Make sure we get dispersion straight from L1B
config.fm.instrument.dispersion.creator = ConfigCommon.dispersion_polynomial
config.fm.instrument.dispersion.apriori = ConfigCommon.l1b_spectral_coefficient_i

config:do_config()
