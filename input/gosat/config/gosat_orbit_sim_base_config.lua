------------------------------------------------------------
--- This configuration modifies the operations configuration
--- for running on orbit simulator data.
---
--- * Change input ILS to have a flat response
--- * Use true noise, not empirical noise
--- * Make ILS width wider
--- * Change prior CO2 by a scale factor
--- * Turn off ZLO
--- * Turn off non-uniform resampling
--- * Turn off aerosols, use rayleigh only (commented out by default)
------------------------------------------------------------

------------------------------------------------------------
--- Start with the norma BaseConfig, this is the normal
--- production run.
------------------------------------------------------------

require "gosat_base_config"

GosatOrbitSimBaseConfig = GosatBaseConfig:new()

---***************************************************************
--- Change input ILS to have a flat response
---***************************************************************

  --- Use an ILS table as the ILS function, but modify the 
  --- response function to be constant across the pixels.

GosatOrbitSimBaseConfig.ils_table_constant_response = Creator:new()

function GosatOrbitSimBaseConfig.ils_table_constant_response:create_i(i)
   --- HDF group to read
   local hdf_group = "Instrument/ILS/ILS_" .. (i + 1) .. "/"
   --- Read data from HDF file
   local wn = self.config:h():read_double_1d(hdf_group .. "wavenumber")
   local del_lam = self.config:h():read_double_2d(hdf_group .. "delta_lambda")
   local response = self.config:h():read_double_2d(hdf_group .. "response")
   local band_name = self.config.common.desc_band_name:value(i)
   local hdf_band_name = self.config.common.hdf_band_name:value(i)
   --- We want to interpolate wavenumbers
   local interp = true

   --- Replace the response with a flat response across all pixels. The 
   --- row we choose depends on i
   for j=0, response:cols() - 1 do
      if(i + 1 == 1) then
	 response:set(0, j, response(1,j))
	 response:set(2, j, response(1,j))
      elseif(i + 1 == 2) then
	 response:set(0, j, response(2,j))
	 response:set(1, j, response(2,j))
	 response:set(3, j, response(2,j))
      elseif(i + 1 == 3) then
	 response:set(1, j, response(0,j))
	 response:set(2, j, response(0,j))
	 response:set(3, j, response(0,j))
      else
	 error("We only support bands 1 through 3")
      end
   end
   return IlsTableLinear(wn, del_lam, response, band_name, hdf_band_name, interp)
end

function GosatOrbitSimBaseConfig.ils_table_constant_response:create()
   local i
   local res = {}
   for i=1,self.config.number_pixel:rows() do
      res[i] = self:create_i(i-1)
   end
   return res
end
GosatOrbitSimBaseConfig.fm.instrument.ils_func.creator = GosatOrbitSimBaseConfig.ils_table_constant_response

---***************************************************************
--- Change to true noise
---***************************************************************

GosatOrbitSimBaseConfig.fm.l1b.noise.creator = AcosConfig.gosat_noise_l1b

---***************************************************************
--- Change ILS width 
---***************************************************************

GosatOrbitSimBaseConfig.fm.instrument.ils_half_width = { DoubleWithUnit(49.99, "cm^-1"),
                                                DoubleWithUnit(49.99, "cm^-1"),
                                                DoubleWithUnit(49.99, "cm^-1") }

---***************************************************************
--- Change prior CO2 by a scale factor
---***************************************************************

   --- We replace the apriori function with this new one function 
   --- with one that scales the original co2 apriori. 

old_aprior_func = GosatOrbitSimBaseConfig.fm.atmosphere.absorber.CO2.apriori
function GosatOrbitSimBaseConfig:scaled_apriori()
   local v = old_aprior_func(self)
   for i=0,v:rows() - 1 do
      v:set(i, v(i) * 0.9873)
   end
   return v
end
GosatOrbitSimBaseConfig.fm.atmosphere.absorber.CO2.apriori = GosatOrbitSimBaseConfig.scaled_apriori

---***************************************************************
--- Change to ABSCO 3.3
---***************************************************************
GosatOrbitSimBaseConfig.fm.atmosphere.absorber.CO2.absco = "v3.2.0/co2_4700.0-6500.0_v3.2.0.chunk.hdf"
GosatOrbitSimBaseConfig.fm.atmosphere.absorber.H2O.absco = "v3.2.0/h2o_4700.0-6500.0_v3.2.0.chunk.hdf"
GosatOrbitSimBaseConfig.fm.atmosphere.absorber.O2.absco = "v3.2.0/o2_12745.0-13245.0_v3.2.0.chunk.hdf"

---***************************************************************
--- Turn off ZLO
---***************************************************************

GosatOrbitSimBaseConfig.fm.instrument.instrument_correction.ic = { }

---***************************************************************
--- Turn off non-uniform sampling
---***************************************************************

GosatOrbitSimBaseConfig.fm.spec_samp.creator = ConfigCommon.uniform_spectrum_sampling

---***************************************************************
--- Turn off aerosols (optionally)
---***************************************************************

--- GosatOrbitSimBaseConfig.fm.atmosphere.aerosol.creator = ConfigCommon.rayleigh_only
