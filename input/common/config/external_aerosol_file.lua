-- Modifications from standard configuration for using a file specifying aerosol apriori information
-- 
-- Requires an external aerosol file formatted like so:
-- /Aerosol/GaussianParam   Dataset {F, S, T, A}
-- /Aerosol/GaussianParamCov Dataset {F, S, T, A, A}
-- /Aerosol/TypeNames       Dataset {T}
-- /Aerosol/TypesUsed       Dataset {F, S, T}
--
-- F = frame index
-- S = sounding index
-- T = aerosol type index
-- A = aerosol param index
--
-- Use by putting the following in a config.lua file:
--
-- config.aer_file = "/path/to/external_aerosol_file.h5"
--
-- require "external_aerosol_file"
-- init_external_aerosol_file(config)
-- 

-- For table.contains
require "helper_functions"

function init_external_aerosol_file(config)

    -------------
    -- Aerosol --
    -------------

    -- For now pass as an argument but move to config later
    if (not config.aer_file) then
        error("aer_file not defined in config")
    end

    -- Define function to read aerosol priori file into the config
    function config:aer_priori_file()
       if(self.aer_file and not self.aer_file_v) then
          self.aer_file_v = HdfFile(self.aer_file)
          self.input_file_description = self.input_file_description .. 
            "Aerosol priori input file:      " .. self.aer_file .. "\n"
       end
       return self.aer_file_v
    end

    -- Overall aerosol creator that sets up the structure for those aerosol types present in the external file
    ext_aerosol_creator = ConfigCommon.aerosol_creator:new()

    function ext_aerosol_creator:sub_object_key()
        local sid = self.config:l1b_sid_list()
        local frame_idx = sid:frame_number()
        local sounding_idx = self.config.sid:sounding_number()
        local type_indexes = self.config:aer_priori_file():read_int_3d("/Aerosol/TypesUsed")(frame_idx, sounding_idx, Range.all()) 
        local type_names = self.config:aer_priori_file():read_string_vector("/Aerosol/TypeNames")
        
        local aerosols = {}
        for idx=0,type_indexes:rows()-1 do
            if (type_indexes(idx) >= 0) then
                local type_idx = type_indexes(idx)
                table.insert(aerosols, type_names:value(type_idx))
            end
        end

        return aerosols
    end

    function config:aerosol_type_index(aer_name)
        local sid = self:l1b_sid_list()
        local frame_idx = sid:frame_number()
        local sounding_idx = sid:sounding_number()
        
        local type_indexes = self:aer_priori_file():read_int_3d("/Aerosol/TypesUsed")(frame_idx, sounding_idx, Range.all()) 
        local type_names = self:aer_priori_file():read_string_vector("/Aerosol/TypeNames")

        local aer_index
        for idx=0,type_indexes:rows()-1 do
            if (type_indexes(idx) >= 0) then
                local type_idx = type_indexes(idx)
                if (type_names:value(type_idx) == aer_name) then
                    aer_index = idx
                end
            end
        end

        if(aer_index == nil) then
            error("Could not find profile index for aerosol type: " .. aer_name)
        end

        return aer_index
    end

    function ext_aerosol_apriori()
        return function(self, aer_name)
            local sid = self.config:l1b_sid_list()
            local frame_idx = sid:frame_number()
            local sounding_idx = sid:sounding_number()

            local type_index = self.config:aerosol_type_index(aer_name)
            return self.config:aer_priori_file():read_double_4d("/Aerosol/GaussianParam")(frame_idx, sounding_idx, type_index, Range.all())
        end
    end

    function ext_aerosol_covariance()
        return function(self, aer_name)
            local sid = self.config:l1b_sid_list()
            local frame_idx = sid:frame_number()
            local sounding_idx = sid:sounding_number()

            local type_index = self.config:aerosol_type_index(aer_name)
            val = self.config:aer_priori_file():read_double_5d("/Aerosol/GaussianParamCov")(frame_idx, sounding_idx, type_index, Range.all(), Range.all())
            print("cov = ", val)
            return val
        end
    end

    config.fm.atmosphere.aerosol = {
        creator = ext_aerosol_creator
    }

    aerosol_prop_names = config:aer_priori_file():read_string_vector("/Aerosol/TypeNames")

    for i = 0, aerosol_prop_names:size()-1 do
        prop_name = aerosol_prop_names:value(i)
        config.fm.atmosphere.aerosol[prop_name] = {
            creator = ConfigCommon.aerosol_log_shape_gaussian,
            apriori = ext_aerosol_apriori(),
            covariance = ext_aerosol_covariance(),
            property = ConfigCommon.hdf_aerosol_property(prop_name),
        }
    end

end
