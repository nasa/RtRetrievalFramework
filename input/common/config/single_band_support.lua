------------------------------------------------------------
--- Enables single band selection for instruments with
--- OCO like band names: ABO2, WCO2, SCO2
------------------------------------------------------------

function init_single_band_support(config)
    ------------------------------------------------------------
    --- Selects which spectrometers from the static HDF file to
    --- use based on the setting of the spectrometer's variable
    --- from the config file. The default is to use all 
    --- spectrometers.
    ------------------------------------------------------------
    orig_spec_win_creator = config.fm.spec_win.creator 
    config.spec_select_spectral_window_hdf = Creator:new()
    function config.spec_select_spectral_window_hdf:create()
        local win_ranges = orig_spec_win_creator.create(self):range_array()
        -- Determine which spectrometers to disable if the which_spectrometers 
        -- variable is defined
        local which_spec = self.config.which_spectrometers
        if which_spec then
            if not string.find(which_spec, "ABO2") then
                win_ranges.value:set(0, 0, Range.all(), 0.0)
            end
            if not string.find(which_spec, "WCO2") then
                win_ranges.value:set(1, 0, Range.all(), 0.0)
            end
            if not string.find(which_spec, "SCO2") then
                win_ranges.value:set(2, 0, Range.all(), 0.0)
            end
        end

        if (self.bad_sample_mask) then
            local bsamp = self:bad_sample_mask()
            if(bsamp) then
                return SpectralWindowRange(win_ranges, bsamp)
            else
                return SpectralWindowRange(win_ranges)
            end
        else
            return SpectralWindowRange(win_ranges)
        end
    end
    config.fm.spec_win.creator = config.spec_select_spectral_window_hdf

    ------------------------------------------------------------
    -- Set which spectrometers will be included in the retrieval
    -- If this is not set then all bands are retrieved. If it is
    -- set then the string should contain one or more of the
    -- following strings:
    -- "ABO2", "WCO2", "SCO2"
    -- Without the commas or quotes of course. For example
    -- The default configuration would be equivalent to:
    -- which_spectrometers = "ABO2 WCO2 SCO2"
    ------------------------------------------------------------
    if not config.which_spectrometers or config.which_spectrometers == "" then
        config.which_spectrometers = os.getenv("which_spectrometers")
    end

    -- Set up options needed for single band retrievals
    if (config.which_spectrometers) then
        -- Turn off surface pressure if when the aband is not included
        if (not string.find(config.which_spectrometers, "ABO2")) then
           config.fm.atmosphere.pressure.retrieved = false
        end

        -- Remove CO2 from gases if not doing CO2 bands
        -- Turn off H2O retrieval
        if ((not string.find(config.which_spectrometers, "WCO2")) and (not string.find(config.which_spectrometers, "SCO2"))) then
           config.fm.atmosphere.absorber.H2O.retrieved = false
           config.fm.atmosphere.absorber.gases = {"H2O", "O2"} 
        end
    end
end
