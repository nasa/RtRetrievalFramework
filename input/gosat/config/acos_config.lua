-- This adds ACOS/GOSAT specific routines for use by base_config

require "config_common"

------------------------------------------------------------
--- Get the file names found in source_files, because
--- they change from run to run.
------------------------------------------------------------

AcosConfig = ConfigCommon:new()

------------------------------------------------------------
-- Determine sounding id list from HDF file
-- Overrides empty behavior from ConfigCommon
------------------------------------------------------------

function AcosConfig:l1b_sid_list()
   if(not self.l1b_sid_list_v) then
      local hv = self:l1b_hdf_file()
      if(hv) then
         self.l1b_sid_list_v = AcosSoundingId.create(hv, self.sid_string)
	 self.sid = self.l1b_sid_list_v:value(0)
      end
   end
   return self.l1b_sid_list_v
end

------------------------------------------------------------
--- Create ECMWF if we have one
------------------------------------------------------------

AcosConfig.acos_met = Creator:new()

function AcosConfig.acos_met:create()
   local sid = self.config:l1b_sid_list()
   if (self.config.met_file) then
       local met = AcosMetFile(self.config.met_file, self.config:l1b_sid_list():value(0),
                               self.config:l1b_sid_list():size() > 1)
       self.config.input_file_description = self.config.input_file_description .. 
          "Meteorology input file:    " .. self.config.met_file .. "\n"
       return met
   end
end

function AcosConfig.acos_met:register_output(ro)
    if (self.config.met) then
        ro:push_back(MetPassThroughOutput(self.config.met))
    end
end

------------------------------------------------------------
--- Create a level 1b file were we read it from a HDF file.
---
---  Depending on the value of sid_string, we either do S, P,
---  or an average of both. So "20090726171701S" uses S band,
---  "20090726171701P" uses P band and "20090726171701" averages.
------------------------------------------------------------

AcosConfig.level1b_hdf = CreatorL1b:new()

function AcosConfig.level1b_hdf:create_parent_object()
   local hv = self.config:l1b_hdf_file()
   local slist = self.config:l1b_sid_list()
   local l1blist = VectorLevel1b()
   for i=0,slist:size() - 1 do
      local l1bacos = Level1bAcos(hv, slist:value(i))
      l1bacos.noise_model = self.config.noise[i]
      l1blist:push_back(l1bacos)
      
      if (self.config.land_fraction and 
          self.config.land_fraction ~= l1bacos:land_fraction(0)) then
         error("Land fraction should be the same for all sounding parts")
      else
         self.config.land_fraction = l1bacos:land_fraction(0)
      end
      self.config.is_h_gain = l1bacos:is_h_gain()
   end
   return Level1bAverage(l1blist)
end

------------------------------------------------------------
--- Variation of spectrum_effect_list that has a different
--- set of spectrum effects for M for H gain.
------------------------------------------------------------

AcosConfig.spectrum_effect_list_h_and_m = ConfigCommon.spectrum_effect_list:new()

function AcosConfig.spectrum_effect_list_h_and_m:sub_object_key()
   local res
   if(self.config.is_h_gain) then
      res = self.speceff_h_gain
   else
      res = self.speceff_m_gain
   end
   return res
end

------------------------------------------------------------
--- Variation of instrument_correction_list that has a different
--- set of spectrum effects for M for H gain.
------------------------------------------------------------

AcosConfig.instrument_correction_list_h_and_m = ConfigCommon.instrument_correction_list:new()

function AcosConfig.instrument_correction_list_h_and_m:sub_object_key()
   local res
   if(self.config.is_h_gain) then
      res = self.ic_h_gain
   else
      res = self.ic_m_gain
   end
   return res
end

------------------------------------------------------------
--- Create a Gosat Noise model.
------------------------------------------------------------

AcosConfig.gosat_noise = Creator:new()

function AcosConfig.gosat_noise:create(hdf_file, sid)
   local hv = self.config:l1b_hdf_file()
   local slist = self.config:l1b_sid_list()
   local res = {}
   for i=0,slist:size() - 1 do
      res[i] = GosatNoiseModel(hv, slist:value(i), self.config.common.hdf_band_name, 
                               self.config:h(), "Instrument")
   end
   return res
end

------------------------------------------------------------
--- Use a Gosat Noise model, but no empirical noise (just
--- what is found in level 1b file)
------------------------------------------------------------

AcosConfig.gosat_noise_l1b = Creator:new()

function AcosConfig.gosat_noise_l1b:create(hdf_file, sid)
   local hv = self.config:l1b_hdf_file()
   local slist = self.config:l1b_sid_list()
   local res = {}
   for i=0,slist:size() - 1 do
      res[i] = GosatNoiseModel(hv, slist:value(i), self.config.common.hdf_band_name)
   end
   return res
end

------------------------------------------------------------
--- Retrieve offset scaling for static HDF file
------------------------------------------------------------

function AcosConfig.gosat_offset_scaling(field)
   return 
      function (self)
         return self.config:h():read_double_with_unit_1d(field)
      end
end

------------------------------------------------------------
--- Return the ground type name 
------------------------------------------------------------

AcosConfig.ground_type_name = DispatchCreator:new()

function AcosConfig:ground_type_name()
   if (not self.land_fraction) then
      error("Land fraction not defined yet")
   end
   if(self.land_fraction >= 0.0 and 
      self.land_fraction <= 90.0) then
      return "coxmunk"
   elseif(self.land_fraction > 90.0 and 
          self.land_fraction <= 100.0) then
      return "lambertian"
   else
      error("Invalid land fraction value read from L1B file: " .. self.config.land_fraction)
   end
end
------------------------------------------------------------
--- Create ground based on the surface type
------------------------------------------------------------

AcosConfig.ground_land_fraction = DispatchCreator:new()

function AcosConfig.ground_land_fraction:get_creator()
   if(self.config:ground_type_name() == "coxmunk") then
      if(self.config.using_radiance_scaling ~= nil) then
         return ConfigCommon.ground_coxmunk
      else
         return ConfigCommon.ground_coxmunk_plus_lamb
      end
   elseif(self.config:ground_type_name() == "lambertian") then
      return ConfigCommon.ground_lambertian
   else
      error("Invalid ground type name: " .. self:ground_type_name())
   end
end

------------------------------------------------------------
--- Create ILS table by reading from an L1B hdf file using
--- Lua to do the reading and index massaging.
------------------------------------------------------------
AcosConfig.ils_table_l1b = Creator:new()

function AcosConfig.ils_table_l1b:create()
   -- Load reference to L1B HDF file
   l1b_hdf_file = self.config:l1b_hdf_file()

   -- O2 band in L1B has an ILS table for P and S, 
   -- Select by default the one used in production (S)
   -- Can be set in base_config.lua
   local o2_table_pick
   if self.o2_table_pick == nil then
      o2_table_pick = 1
   else
      o2_table_pick = self.o2_table_pick
   end

   -- Whether or not to interpolate the table values
   -- Can be set in base_config.lua
   local interpolate
   if self.interpolate == nil then
      interpolate = true
   else
      interpolate = self.interpolate
   end

   -- We know this is an ACOS file, since we had to write
   -- specifics for it, so might as well define band names here
   local idx
   local wavenumber, delta_lambda, response_name
   local res = {}
   for idx = 0, self.config.common.hdf_band_name:size() - 1 do
      local hdf_band_name = self.config.common.hdf_band_name:value(idx)
      local desc_band_name = self.config.common.desc_band_name:value(idx)

      -- Load center wavenumbers of ILS
      if hdf_band_name == "o2" then
         wavenumber = l1b_hdf_file:read_double_2d("InstrumentHeader/ils_coef_center_wavenumber_" .. hdf_band_name)
         wavenumber = wavenumber(o2_table_pick, Range())
      else
         wavenumber = l1b_hdf_file:read_double_1d("InstrumentHeader/ils_coef_center_wavenumber_" .. hdf_band_name)
      end

      -- Load ILS response values
      if hdf_band_name == "o2" then
         response = l1b_hdf_file:read_double_3d("InstrumentHeader/ils_coef_" .. hdf_band_name)
         response = response(o2_table_pick, Range(), Range())
      else
         response = l1b_hdf_file:read_double_2d("InstrumentHeader/ils_coef_" .. hdf_band_name)
      end

      -- Delta lambda as provided in L1B file is meant for all center wavenumbers,
      -- But need to copy that array into a 2d array to meet IlsTable interface
      -- which expects a different delta lambda per wavenumber index
      delta_lambda = Blitz_double_array_2d(response:rows(), response:cols())

      delta_lambda_in = l1b_hdf_file:read_double_1d("InstrumentHeader/ils_coef_relative_wavenumber_" .. hdf_band_name)

      for wn_idx = 0, response:rows() - 1 do
         delta_lambda:set(wn_idx, Range(), delta_lambda_in)
      end

      res[idx+1] = IlsTableLinear(wavenumber, delta_lambda, response, desc_band_name, hdf_band_name, interpolate)
   end
   return res
end

------------------------------------------------------------
--- Determine of we are land or water
------------------------------------------------------------

function AcosConfig:land_or_water()
   if(self.land_fraction > 90.0) then
      return "land"
   else
      return "water"
   end
end

------------------------------------------------------------
--- Fluorescence only for land runs
------------------------------------------------------------

AcosConfig.fluorescence_effect_land_only = ConfigCommon.fluorescence_effect:new()
function AcosConfig.fluorescence_effect_land_only:create()
   -- Only use this for land runs
   if(self.config:land_or_water() == "land") then
      return ConfigCommon.fluorescence_effect.create(self)
   else
      return nil
   end
end

------------------------------------------------------------
--- Defines the ground type name we should use, based on
--- land water fag
------------------------------------------------------------

function AcosConfig:ground_type_name()
   if(self:land_or_water() == "land") then
      return "brdf_soil"
   else
      return "coxmunk"
   end
end

------------------------------------------------------------
--- Create ground based on the surface type
------------------------------------------------------------

AcosConfig.ground_from_ground_type = DispatchCreator:new()

function AcosConfig.ground_from_ground_type:get_creator()
   local ground_type = self.config:ground_type_name()

   if (ground_type == "lambertian") then
      return ConfigCommon.ground_lambertian
   elseif (ground_type == "brdf_soil") then
      return ConfigCommon.ground_brdf_soil
   elseif (ground_type == "brdf_veg") then
      return ConfigCommon.ground_brdf_veg
   elseif (ground_type == "coxmunk") then
      if(self.config.using_radiance_scaling ~= nil) then
         return ConfigCommon.ground_coxmunk
      else
         return ConfigCommon.ground_coxmunk_plus_lamb
      end
   else
      error("Invalid ground type value: " .. ground_type)
   end
end


------------------------------------------------------------
--- Use tropopause height for initial guess as cirrus ice
--- height
------------------------------------------------------------

function AcosConfig.tropopause_height_ap(self, base, type, aer_name)
    local ap = self:apriori_initial(base, type, aer_name)
    
    local t = self.config:reference_co2_apriori_met_obj()
    local psurf = self.config.met:surface_pressure()
    ap:set(1, (t:tropopause_pressure() - 100) / psurf)
    
    return ap
end
