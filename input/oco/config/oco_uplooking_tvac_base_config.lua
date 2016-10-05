------------------------------------------------------------
-- Makes the  necessary changes to the nominal OCO
-- configuration to support TVAC2 uplooking data retrievals.
------------------------------------------------------------

require "oco_uplooking_base_config"

tvac_config = OcoUplookingBaseConfig:new()

-- Fix dispersion covariance
orig_disp_cov = tvac_config.fm.instrument.dispersion.covariance
function tvac_config.dispersion_cov_size_6(self, i)
   cov = orig_disp_cov(self, i)
   cov = cov(Range(0,5), Range(0,5))
   return cov
end
tvac_config.fm.instrument.dispersion.covariance = tvac_config.dispersion_cov_size_6

-- Use a surface pressure covariance better suited for tvac runs
function tvac_press_cov(self)
    cov = Blitz_double_array_2d(1,1)
    cov:set(0, 0, 2.5e7)
    return cov
end

tvac_config.fm.atmosphere.pressure.covariance = tvac_press_cov

-- Constrain CO2 profile with a covariance matrix that simulates a scaling retrieval
-- This function is also needed since the CO2 covariance by default is 20 x 20
-- and might not be sized according to the one coming in from a text file
function tvac_co2_cov(self)
    co2_ap = tvac_config.fm.atmosphere.absorber.CO2.apriori(self)
    cov = Blitz_double_array_2d(co2_ap:rows(), co2_ap:rows())
    cov:set(Range.all(), Range.all(), 2.499e-9)
    for idx = 0, cov:rows() -1 do
        cov:set(idx, idx,  2.5e-9)
    end
    return cov
end
tvac_config.fm.atmosphere.absorber.CO2.covariance = tvac_co2_cov

