------------------------------------------------------------
--- Adds to the L2 retrieval additional interference
--- aerosol types that are present in the statevector
--- but have their covariances set such that they don't
--- update. The end result is that these types have
--- their jacobians calculated so that offline 
--- error analysis can be run to determine aerosol 
--- interference.
---
--- This config requires that LIDORT 3.5 be compiled with
--- MAX_ATMOSWFS = 16
------------------------------------------------------------

function init_aerosol_interference(config)

    -- Use a custom setup creator for linear aerosols
    -- to shut off unneeded jacobians for interference aerosols
    config.aer_interf_linear_shape = ConfigCommon.aerosol_linear_shape_gaussian:new()

    function config.aer_interf_linear_shape:retrieval_flag() 
        return ConfigCommon.create_flag(self:apriori_v():rows(), Range(0,0))
    end

    function interference_apriori_1()
        local ap = Blitz_double_array_1d(3)
        -- 0.0, 0.95, 0.05
        ap:set(0, 1.0e-4)
        ap:set(1, 0.95)
        ap:set(2, 0.05)
        return ap
    end

    function interference_apriori_2()
        local ap = Blitz_double_array_1d(3)
        -- 0.0, 0.5, 0.2
        ap:set(0, 1.0e-4)
        ap:set(1, 0.5)
        ap:set(2, 0.2)
        return ap
    end

    function interference_covariance()
        local cov = Blitz_double_array_2d(3,3)
        cov:set(Range.all(), Range.all(), 0.0)
        for idx = 0, 2 do
            cov:set(idx, idx, 1e-20)
        end
        return cov
    end

    -- Add interference aerosols
    local interf_aerosols_propname = { "adu", "ass", "asu", "bc", "oc" } 
    local interf_aerosols_dispname = { "DU", "SS", "SO", "BC", "OC" }
    local interf_ig_groups = { interference_apriori_1, interference_apriori_2 }
    for aer_idx, aer_name in ipairs(interf_aerosols_dispname) do
        local prop_name = interf_aerosols_propname[aer_idx]
        for ig_idx, ig_func in ipairs(interf_ig_groups) do
            local config_aer_name = aer_name .. "_Interference_" .. ig_idx
            table.insert(config.fm.atmosphere.aerosol.aerosols, config_aer_name) 
            config.fm.atmosphere.aerosol[config_aer_name] = {
               creator = config.aer_interf_linear_shape,
               apriori = ig_func,
               covariance = interference_covariance,
               property = ConfigCommon.hdf_aerosol_property(prop_name),
            }
        end
    end

end
