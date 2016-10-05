------------------------------------------------------------
--- Performs a FTS simulation only
------------------------------------------------------------

require "fts_base_config"

config = FtsBaseConfig:new()

-- Diagnostics
config.diagnostic = false

-- Use only WCO2 band
spectral_window_wco2_only = Creator:new()

function spectral_window_wco2_only:create()
   orig_range = self.config:h():read_double_with_unit_3d("Spectral_Window/microwindow")
   -- "zero" out ABO2 and SCO2 bands
   orig_range.value:set(0, 0, Range.all(), 0.0)
   orig_range.value:set(2, 0, Range.all(), 0.0)
   return SpectralWindowRange(orig_range)
end
config.fm.spec_win.creator = spectral_window_wco2_only

-- Do not retrieve surface pressure
config.fm.atmosphere.pressure.retrieved = false

-- Set up file paths
config.do_retrieval = false
config.spectrum_1_file = ""
config.spectrum_2_file = "input/pa20091103saaaaa_100223160344.008"
config.spectrum_3_file = ""
config.runlog_file = "input/tccon_runlog.grl"
config.atmosphere_file = "input/atmosphere_fts.dat"

config:do_config()
 
