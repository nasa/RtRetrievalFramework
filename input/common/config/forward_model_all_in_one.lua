-- Modifications from standard configuration for a forward model input file
-- that contains all necessary processing information

-- For table.contains
require "helper_functions"

function init_forward_model(config)

    config.diagnostic = true

    -- Turn of retrieval mode
    config.do_retrieval = False

    -- End after one iteration of the mechanisms
    config.solver.max_iteration = 1

    -- Output jacobians
    config.write_jacobian = true

    ------------
    -- Common --
    ------------
    
    function sim_apriori_i(field)
        return function(self, i)
            return self.config:l1b_hdf_file():read_double_1d(field .. "/a_priori")(i)
        end
    end

    function sim_apriori_sounding(field)
        return function(self, aer_name)
            local sid = self.config:l1b_sid_list()
            local frame_idx = sid:frame_number()
            local sounding_idx = sid:sounding_number()

            return self.config:l1b_hdf_file():read_double_3d(field)(frame_idx, sounding_idx, Range.all())
        end
    end

    ---------------------
    -- Spectral Window --
    ---------------------
    
    local orig_spec_win_creator = config.fm.spec_win.creator

    spectral_window_custom = Creator:new()
    function spectral_window_custom:create()
        local orig_win_ranges = orig_spec_win_creator.create(self):range_array()

        local new_win_data = Blitz_double_array_3d(orig_win_ranges.value:rows(), 1, 2)
        new_win_data:set(Range.all(), Range.all(), Range.all(), 0.0)
        local new_win_ranges = ArrayWithUnit_3d(new_win_data, orig_win_ranges.units)

        for band_idx = 0, new_win_data:rows()-1 do
            new_win_ranges.value:set(band_idx, 0, 0, 1)
            new_win_ranges.value:set(band_idx, 0, 1, 1017)
        end

        if (self.bad_sample_mask) then
            local bsamp = self:bad_sample_mask()
            if(bsamp) then
                return SpectralWindowRange(new_win_ranges, bsamp)
            else
                return SpectralWindowRange(new_win_ranges)
            end
        else
            return SpectralWindowRange(new_win_ranges)
        end

    end
    config.fm.spec_win.creator = spectral_window_custom

    -- Use bad samples mask only, don't use SAA masking
    config.fm.spec_win.bad_sample_mask = OcoConfig.snr_coef_bad_sample_mask

    --------------
    -- Pressure --
    --------------

    config.fm.atmosphere.pressure.creator = ConfigCommon.pressure_sigma_profile
    config.fm.atmosphere.pressure.pressure_levels = sim_apriori_sounding("/Pressure/PressureLevels")

    ------------
    -- Ground --
    ------------
    config.fm.atmosphere.ground.lambertian.apriori = sim_apriori_i("Ground/Albedo")

    sim_ground_selector = DispatchCreator:new()

    function sim_ground_selector:get_creator()
       local ground_type = self.config:ground_type_name()

       if (ground_type == "lambertian") then
          return ConfigCommon.ground_lambertian
       elseif (ground_type == "brdf_soil") then
          return ConfigCommon.ground_brdf_soil
       elseif (ground_type == "brdf_veg") then
          return ConfigCommon.ground_brdf_veg
       elseif (ground_type == "coxmunk") then
          return ConfigCommon.ground_coxmunk
       else
          error("Invalid ground type value: " .. ground_type)
       end
    end

    config.fm.atmosphere.ground.creator = sim_ground_selector

    -------------
    -- Aerosol --
    -------------
    sim_aerosol_creator = ConfigCommon.aerosol_creator:new()

    function sim_aerosol_creator:sub_object_key()
        local sid = self.config:l1b_sid_list()
        local frame_idx = sid:frame_number()
        local sounding_idx = self.config.sid:sounding_number()
        local type_indexes = self.config:l1b_hdf_file():read_int_3d("/Aerosol/TypesUsed")(frame_idx, sounding_idx, Range.all()) 
        local type_names = self.config:l1b_hdf_file():read_string_vector("/Aerosol/TypeNames")
        
        local aerosols = {}
        for idx=0,type_indexes:rows()-1 do
            if (type_indexes(idx) >= 0) then
                type_idx = type_indexes(idx)
                table.insert(aerosols, type_names:value(type_idx))
            end
        end

        return aerosols
    end

    sim_aerosol_profile_creator = ConfigCommon.aerosol_linear_profile:new()

    function sim_aerosol_profile_creator:covariance_v()
        local cov = Blitz_double_array_2d(20, 20)
        cov:set(Range.all(), Range.all(), 1e-20)
        return cov
    end

    function sim_aerosol_apriori()
        return function(self, aer_name)
            local sid = self.config:l1b_sid_list()
            local frame_idx = sid:frame_number()
            local sounding_idx = sid:sounding_number()
            
            local type_indexes = self.config:l1b_hdf_file():read_int_3d("/Aerosol/TypesUsed")(frame_idx, sounding_idx, Range.all()) 
            local type_names = self.config:l1b_hdf_file():read_string_vector("/Aerosol/TypeNames")

            for idx=0,type_indexes:rows()-1 do
                if (type_indexes(idx) >= 0) then
                    type_idx = type_indexes(idx)
                    if (type_names:value(type_idx) == aer_name) then
                        profile_idx = idx
                    end
                end
            end
            if(profile_idx == nil) then
                error("Could not find profile index for aerosol type: " .. aer_name)
            end

            return self.config:l1b_hdf_file():read_double_4d("/Aerosol/Profiles")(frame_idx, sounding_idx, profile_idx, Range.all())
        end
    end

    config.fm.atmosphere.aerosol = {
        creator = sim_aerosol_creator
    }

    aerosol_prop_types ={ "wc_004", "wc_005", "wc_006", "wc_007", "wc_008", "wc_009", "wc_010",
        "wc_011", "wc_012", "wc_013", "wc_014", "wc_015", "wc_016", "wc_017",
        "wc_018", "wc_019", "wc_020", "wc_021", "wc_022", "wc_023",
        "wc_024", "wc_025", "wc_026", "wc_027", "wc_028", "wc_029", "wc_030",
        "ic_010", "ic_015", "ic_020", "ic_025", "ic_030", "ic_035", "ic_040",
        "ic_045", "ic_050", "ic_055", "ic_060", "ic_065", "ic_070", "ic_075",
        "ic_080", "ic_085", "ic_090",
    }

    for i, prop_name in ipairs(aerosol_prop_types) do
        config.fm.atmosphere.aerosol[prop_name] = {
            creator = sim_aerosol_profile_creator,
            apriori = sim_aerosol_apriori(),
            property = ConfigCommon.hdf_aerosol_property(prop_name),
        }
    end

    -- Remove EOF from state vector
    config.fm.instrument.instrument_correction.ic_nadir = {}
    config.fm.instrument.instrument_correction.ic_glint = {}
    config.fm.instrument.instrument_correction.ic_target = {}

end
