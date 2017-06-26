------------------------------------------------------------
--- Runs the BRDF case as if it were a lambertian using
--- the BRDF kernel
------------------------------------------------------------

require "oco_base_config"

config = OcoBaseConfig:new()

-- Use BRDF Soil Ground
use_brdf = DispatchCreator:new()

function use_brdf:get_creator()
    return ConfigCommon.ground_brdf_soil
end

config.fm.atmosphere.ground.creator = use_brdf

function brdf_retrieval_flag(self, i)
    local flag = Blitz_bool_array_1d(self:apriori_v(i - 1):rows())

    flag:set(Range.all(), false)
    -- Weight intercept
    flag:set(0, true)
    -- Weight slope
    flag:set(1, true)

    return flag
end

config.fm.atmosphere.ground.brdf_veg.retrieval_flag = brdf_retrieval_flag
config.fm.atmosphere.ground.brdf_soil.retrieval_flag = brdf_retrieval_flag

function lambertian_like_ap(self, i)
    lamb_ap = OcoConfig.oco_albedo_from_radiance(0)(self, i)

    local brdf_ap = Blitz_double_array_1d(7)

    brdf_ap:set(Range.all(), 0)

    -- Overall weight intercept
    brdf_ap:set(0, lamb_ap(0))
    -- Rahman factor
    brdf_ap:set(2, 1.0)
    -- Overall amplitude
    brdf_ap:set(3, 1.0)
    -- Asymmetry parameter
    brdf_ap:set(4, 0.0)
     -- Geometric factor
    brdf_ap:set(5, 1.0)
    -- Breon factor
    brdf_ap:set(6, 1.0e-20)

    return brdf_ap
end

config.fm.atmosphere.ground.brdf_soil.apriori = lambertian_like_ap
config.fm.atmosphere.ground.brdf_veg.apriori = lambertian_like_ap

config:do_config()
