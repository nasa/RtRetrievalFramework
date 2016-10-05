------------------------------------------------------------
--- This is the configuration used in the "quick look" test
--- (the 25,000 soundings test). 
---
--- If these changes get approved, they will need to be moved
--- into the production configuration file (e.g., gosat_base_config),
--- but we keep these separate while testing
------------------------------------------------------------

require "gosat_base_config"

config = GosatBaseConfig:new()

---------------------------------------------------------------
--- Qlert #0, run with 150/151 data. We have updated the EOFs 
--- for 160/161, but have the old EOF still available
---------------------------------------------------------------

-- Note the apriori and covariance is unchanged between the old and new
-- EOF, so we don't need to bother changing them. Only the hdf group used
-- to read the actual EOF needs to be changed.
-- config.fm.instrument.instrument_correction.eof_h_gain_1.hdf_group = "Instrument/Eof_150151/H_gain"
-- config.fm.instrument.instrument_correction.eof_m_gain_1.hdf_group = "Instrument/Eof_150151/M_gain"

---------------------------------------------------------------
--- Qlert #0.5, use 160/161 data. No changes, just comment out
--- the use of the 150151 EOFs above
---------------------------------------------------------------

---------------------------------------------------------------
--- Qlert #1, turn off EOFs
---------------------------------------------------------------

-- config.fm.instrument.instrument_correction.ic_h_gain = {"zero_offset_waveform"}
-- config.fm.instrument.instrument_correction.ic_m_gain = {"zero_offset_waveform"}

---------------------------------------------------------------
--- Qlert #2 Add methane
---------------------------------------------------------------

-- config.fm.atmosphere.absorber.gases = {"CO2", "H2O", "O2", "CH4"}
-- config.fm.atmosphere.absorber.CH4 = {
--    apriori = function(self)
--       local vmr = Blitz_double_array_1d(20)
--       for i=1,20 do
-- 	 vmr:set(i-1, 1.8e-6)
--       end
--       return vmr
--    end,
--    absco = "v4.3.0_testing_unstable/narrow/ch4_hitran2012.hdf",
--    creator = ConfigCommon.vmr_level_constant_well_mixed,
-- }

---------------------------------------------------------------
--- Qlert #3 New H2O list
---------------------------------------------------------------

-- config.fm.atmosphere.absorber.H2O.absco = "v4.3.0_testing_unstable/narrow/h2o_hitran2012.hdf"

---------------------------------------------------------------
--- Qlert #4 Update to the 1.6 micron CO2 band
---------------------------------------------------------------

-- config.fm.atmosphere.absorber.CO2.absco = "v4.3.0_testing_unstable/narrow/co2_v4.2.0_updated_wco2.hdf"
-- config.fm.atmosphere.absorber.CO2.table_scale = {1.0, 1.0078, 0.9946}

---------------------------------------------------------------
--- Qlert #5 Update to the 2.06 micron SCO2 band
---------------------------------------------------------------

-- config.fm.atmosphere.absorber.CO2.absco = "v4.3.0_testing_unstable/narrow/co2_v4.2.0_updated_wco2_updated_sco2.hdf"
-- config.fm.atmosphere.absorber.CO2.table_scale = {1.0, 1.0078, 1.0006}

---------------------------------------------------------------
--- Qlert #6 Remove empirical “CO2 CIA” from the SCO2 band
---------------------------------------------------------------

-- config.fm.atmosphere.absorber.CO2.absco = "v4.3.0_testing_unstable/narrow/co2_v4.2.0_updated_wco2_updated_sco2_noSco2Ctm.hdf"
-- config.fm.atmosphere.absorber.CO2.table_scale = {1.0, 1.0078, 1.0016}

---------------------------------------------------------------
--- Qlert #7 Fit a polynomial albedo, P-branch only for SCO2 band
---------------------------------------------------------------

--- Have 3rd order polynomial fit, but only for strong band.
-- config.fm.atmosphere.ground.lambertian.apriori =
--    function (self, i)
--       local porder = 1
--       -- For SCO2 band, use order 3 polynomial
--       if(i == 2) then
-- 	 porder = 3
--       end
--       return AcosConfig.gosat_albedo_from_radiance(porder)(self, i)
--    end

--- Have normal spectral window, except for SCO2 that is P-branch only

-- spectral_window_p_branch = Creator:new()

-- function spectral_window_p_branch:create()
--     local win_ranges = self.config:h():read_double_with_unit_3d("Spectral_Window/microwindow")

--     win_ranges.value:set(2,0,1, 4850.0)
--     if (self.bad_sample_mask and self:bad_sample_mask()) then
--         return SpectralWindowRange(win_ranges, self:bad_sample_mask())
--     else
--         return SpectralWindowRange(win_ranges)
--     end
-- end

-- config.fm.spec_win.creator = spectral_window_p_branch

---------------------------------------------------------------
--- Qlert #8 Fit a polynomial albedo, R-branch only for SCO2 band
---------------------------------------------------------------

--- Have normal spectral window, except for SCO2 that is R-branch only

-- spectral_window_r_branch = Creator:new()

-- function spectral_window_r_branch:create()
--     local win_ranges = self.config:h():read_double_with_unit_3d("Spectral_Window/microwindow")

--     win_ranges.value:set(2,0,0, 4850.0)
--     win_ranges.value:set(2,0,1, 4888.0)
--     if (self.bad_sample_mask and self:bad_sample_mask()) then
--         return SpectralWindowRange(win_ranges, self:bad_sample_mask())
--     else
--         return SpectralWindowRange(win_ranges)
--     end
-- end

-- config.fm.spec_win.creator = spectral_window_r_branch

---------------------------------------------------------------
--- Qlert #9 EOF updated for QLERT data
---------------------------------------------------------------

-- Note the apriori and covariance is unchanged between the old and new
-- EOF, so we don't need to bother changing them. Only the hdf group used
-- to read the actual EOF needs to be changed.

config.fm.instrument.instrument_correction.eof_h_gain_1.hdf_group = "Instrument/EmpiricalOrthogonalFunctionTest/H_gain"
config.fm.instrument.instrument_correction.eof_m_gain_1.hdf_group = "Instrument/EmpiricalOrthogonalFunctionTest/M_gain"

config:do_config()


