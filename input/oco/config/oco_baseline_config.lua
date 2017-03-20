------------------------------------------------------------
--- This is the configuration used to run against orbit
--- simulator data. 
---
--- There are a number of changes for the normal production run:
--- Disable the ones that are unneeded.
--- * Change prior CO2 by a scale factor
--- * Change ABSCO verison
------------------------------------------------------------

require "oco_base_config"

OcoBaselineConfig = OcoBaseConfig:new()

---***************************************************************
--- Change prior CO2 by a scale factor
---***************************************************************

--- We replace the apriori function with this new one function 
--- with one that scales the original co2 apriori. 
--- Scale factor is determined by comparing the mean true XCO2
--- from Orbit Simulator log files with the L2 calculated xco2_apriori
--- for several thousand soundings. The ratio is what it takes to scale
--- L2's xco2_apriori mean to match that of the log file's XCO2 value

old_aprior_func = OcoBaselineConfig.fm.atmosphere.absorber.CO2.apriori
function OcoBaselineConfig:scaled_apriori()
   local v = old_aprior_func(self)
   for i=0,v:rows() - 1 do
      v:set(i, v(i) * 0.9812)
   end
   return v
end
OcoBaselineConfig.fm.atmosphere.absorber.CO2.apriori = OcoBaselineConfig.scaled_apriori

-- Use CSU ECMWF format reader
OcoBaselineConfig.fm.input.met.creator = OcoConfig.oco_meteorology
