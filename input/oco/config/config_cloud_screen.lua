------------------------------------------------------------
--- Configuration sort of like how the CWO cloud screener
--- works.
------------------------------------------------------------

require "oco_baseline_config"

config = OcoBaselineConfig:new()

config.diagnostic = true

-- Use ABO2 only
spectral_window_abo2_only = Creator:new()
function spectral_window_abo2_only:create()
   orig_range = self.config:h():read_double_with_unit_3d("Spectral_Window/microwindow")
   -- "zero" out WCO2 and SCO2 bands
   orig_range.data:set(1, 0, Range.all(), 0.0)
   orig_range.data:set(2, 0, Range.all(), 0.0)
   return SpectralWindowRange(orig_range)
end
config.fm.spec_win.creator = spectral_window_abo2_only

-- Use lambertian no mater what ground type
ground_always_lambertian = DispatchCreator:new()
function ground_always_lambertian:get_creator()
   return ConfigCommon.ground_lambertian
end
config.fm.atmosphere.ground.creator = ground_always_lambertian

-- Make sure lambertian albedo never > 1.0, since
-- it could be since we are throwing all scenes at this
-- configuration
orig_lambertian_gsd = OcoConfig.oco_ground_from_radiance("Ground/Lambertian")
function lambertian_upper_limit_gsd(self)
   gsd = orig_lambertian_gsd(self)
   for spec = 0, gsd:apriori():rows() - 1 do
      if gsd:apriori()(spec, 0, 0) > 0.90 then
	 gsd:apriori():set(spec, 0, 0, 0.90)
	 gsd:initial_guess():set(spec, 0, 0, 0.90)
      end
   end
   return gsd
end
config.fm.atmosphere.ground.lambertian_gsd = lambertian_upper_limit_gsd

-- No CO2 since we only have ABO2 band
config.fm.atmosphere.absorber.gases = { "H2O", "O2" }

-- No H2O scaling
config.fm.atmosphere.absorber.H2O.retrieved = false

-- No aerosols
config.fm.atmosphere.aerosol.creator = ConfigCommon.rayleigh_only

-- No radiance scaling
config.fm.instrument.instrument_correction.ic = { }

-- Turn on dispersion retrieval
config.fm.instrument.dispersion.retrieved = true

-- Speed up the RT
--config.fm.rt.lsi_constant.high_stream = 2

config:do_config()
