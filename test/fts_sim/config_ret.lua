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
config.do_retrieval = true
config.spectrum_1_file = "out_sim.expected.h5"
config.spectrum_2_file = ""
config.spectrum_3_file = ""
config.runlog_file = ""
config.atmosphere_file = "input/atmosphere_fts.dat"

config.fm.input.l1b.creator = FtsConfig.level1b_fts_hdf

config:do_config()
 
