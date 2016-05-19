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

        return SpectralWindowRange(new_win_ranges)
    end
    config.fm.spec_win.creator = spectral_window_custom

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
        creator = sim_aerosol_creator,
        wc_004 = {
            creator = sim_aerosol_profile_creator,
            apriori = sim_aerosol_apriori(),
            property = ConfigCommon.hdf_aerosol_property("wc_004"),
        },
        wc_008 = {
            creator = sim_aerosol_profile_creator,
            apriori = sim_aerosol_apriori(),
            property = ConfigCommon.hdf_aerosol_property("wc_008"),
        },
        wc_012 = {
            creator = sim_aerosol_profile_creator,
            apriori = sim_aerosol_apriori(),
            property = ConfigCommon.hdf_aerosol_property("wc_012"),
        },
        wc_016 = {
            creator = sim_aerosol_profile_creator,
            apriori = sim_aerosol_apriori(),
            property = ConfigCommon.hdf_aerosol_property("wc_016"),
        },
        wc_018 = {
            creator = sim_aerosol_profile_creator,
            apriori = sim_aerosol_apriori(),
            property = ConfigCommon.hdf_aerosol_property("wc_018"),
        },
        wc_022 = {
            creator = sim_aerosol_profile_creator,
            apriori = sim_aerosol_apriori(),
            property = ConfigCommon.hdf_aerosol_property("wc_022"),
        },
        wc_026 = {
            creator = sim_aerosol_profile_creator,
            apriori = sim_aerosol_apriori(),
            property = ConfigCommon.hdf_aerosol_property("wc_026"),
        },
        wc_030 = {
            creator = sim_aerosol_profile_creator,
            apriori = sim_aerosol_apriori(),
            property = ConfigCommon.hdf_aerosol_property("wc_030"),
        },
     }

    -- Remove EOF from state vector
    config.fm.instrument.instrument_correction.ic_nadir = {}
    config.fm.instrument.instrument_correction.ic_glint = {}
    config.fm.instrument.instrument_correction.ic_target = {}

end
