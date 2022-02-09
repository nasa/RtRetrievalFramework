
-- This adds OCO specific routines for use by base_config

require "config_common"

------------------------------------------------------------
--- Get the file names found in source_files, because
--- they change from run to run.
------------------------------------------------------------

OcoConfig = ConfigCommon:new()

------------------------------------------------------------
--- Determine acquisition mode (nadir, glint, or target), and
--- select the instrument corrections based on this.
------------------------------------------------------------


OcoConfig.instrument_correction_list_acquisition_mode = ConfigCommon.instrument_correction_list:new()


function OcoConfig.instrument_correction_list_acquisition_mode:sub_object_key()

   local acq_mode 
   if (self.config.unscaled_l1b ~= nil) then
      acq_mode = self.config.unscaled_l1b:acquisition_mode()
   else
      acq_mode = self.config.l1b:acquisition_mode()
   end

   local res
   if acq_mode == "Glint" then 
      res = self.ic_glint
   elseif acq_mode == "Target" then 
      res = self.ic_target
   elseif acq_mode == "Transition" then 
      --- We treat this the same as target. We can make a separate ic_transition
      --- if needed in the future, but for now just do what we do with target
      --- mode
      res = self.ic_target
   elseif acq_mode == "Nadir" then 
      res = self.ic_nadir
   elseif acq_mode == "Science" then
      --- This is the OCO-3 science mode, treated the same as target
      res = self.ic_target
   else
      error("Unrecognized acquistion mode " .. acq_mode)
   end
   return res
end

------------------------------------------------------------
--- Determine land vs water, and
--- select the instrument corrections based on this.
------------------------------------------------------------


OcoConfig.instrument_correction_list_land_water = ConfigCommon.instrument_correction_list:new()


function OcoConfig.instrument_correction_list_land_water:sub_object_key()

   local lw
   lw = self.config:land_or_water()
   
   local res
   if lw == "land" then 
      res = self.ic_land
   elseif lw == "water" then 
      res = self.ic_water
   else
      error("Unrecognized land/water " .. lw)
   end
   return res
end

------------------------------------------------------------
-- Determine sounding id list from HDF file
-- Overrides empty behavior from ConfigCommon
------------------------------------------------------------

function OcoConfig:l1b_sid_list()
    if (not self.sid_string or self.sid_string == "") then
        error("Sounding id string not defined")
    end
   if(not self.l1b_sid_list_v) then
      local hv = self:l1b_hdf_file()
      if(hv) then
         self.l1b_sid_list_v = OcoSoundingId(hv, self.sid_string)
	 self.sid = self.l1b_sid_list_v
      end
   end
   return self.l1b_sid_list_v
end

------------------------------------------------------------
--- Create the ECMWF  fileif we have one
------------------------------------------------------------

OcoConfig.oco_met = Creator:new()

function OcoConfig.oco_met:create()
   local sid = self.config:l1b_sid_list()
   if (self.config.met_file) then
       local met = OcoMetFile(self.config.met_file, self.config:l1b_sid_list())
       self.config.input_file_description = self.config.input_file_description .. 
          "Meteorology input file:    " .. self.config.met_file .. "\n"
       return met
   end
end

function OcoConfig.oco_met:register_output(ro)
    if (self.config.met) then
        ro:push_back(MetPassThroughOutput(self.config.met))
    end
end

------------------------------------------------------------
--- Create a level 1b file were we read it from a HDF file.
------------------------------------------------------------

OcoConfig.level1b_hdf = CreatorL1b:new()

function OcoConfig.level1b_hdf:create_parent_object()
   local hv = self.config:l1b_hdf_file()
   local sid = self.config:l1b_sid_list()
   l1b_oco = Level1bOco(hv, sid)
   l1b_oco.noise_model = self.config.noise
   return l1b_oco
end

------------------------------------------------------------
--- Short routine to mark all pixels as bad if they are outside of the
--- normal spectral window.
---
--- This can be used in combination with bad_sample_noise_model
--- to get the full spectrum.
------------------------------------------------------------

function bad_sample_mask_outside_window(self)
   local bad_sample_mask = self:bad_sample_mask_before_window()
   local window = {{90,925}, {211, 863}, {100, 911}}
   for i=1,3 do
      local lb = window[i][1]
      local ub = window[i][2]
      for j=0,bad_sample_mask:cols() - 1 do
         if(j < lb or j >= ub) then
            bad_sample_mask:set(i - 1, j, 1)
         end
      end
   end
   return bad_sample_mask
end

------------------------------------------------------------
--- Create a noise model where we add a large uncertainty
--- where we have bad samples. This is an alternative to just
--- remove the bad samples from the calculation, useful when we
--- want the residuals of the bad samples calculated but not
--- used in the retrieval (e.g., we are analyzing ARPs).
------------------------------------------------------------

OcoConfig.bad_sample_noise_model = Creator:new()

function OcoConfig.bad_sample_noise_model:create()
   local noise_model = self.creator_before_bad_sample.create(self)
   local bad_sample = self:bad_sample_mask()
   return BadSampleNoiseModel(noise_model, bad_sample, 
			      self.bad_sample_uncertainty)
end

------------------------------------------------------------
--- Create a OCO Noise model.
------------------------------------------------------------

OcoConfig.oco_noise = Creator:new()

function OcoConfig.oco_noise:create()
   local hv = self.config:l1b_hdf_file()
   local sid = self.config:l1b_sid_list()
   local nspec = self.config.spec_win:number_spectrometer()
   local max_ms = Blitz_double_array_1d(nspec)
   for i, v in ipairs(self.max_ms) do
      max_ms:set(i-1, v)
   end
   return OcoNoiseModel(hv, sid, max_ms)
end

------------------------------------------------------------
--- Get bad sample from either snr_coef or bad_sample_list,
--- depending on what we find in the file.
------------------------------------------------------------

function OcoConfig.l1b_bad_sample_mask(self)
   local l1b_hdf_file = self.config:l1b_hdf_file()
   if(l1b_hdf_file:has_object("/InstrumentHeader/bad_sample_list")) then
      return self.config.bad_sample_list_bad_sample_mask(self)
   else
      return self.config.snr_coef_bad_sample_mask(self)
   end
end   

------------------------------------------------------------
--- Retrieve a bad pixel mask from the bad_sample_list field.
--- This was added in B8.00
------------------------------------------------------------

function OcoConfig.bad_sample_list_bad_sample_mask(self)
    self.config:diagnostic_message("Using bad_sample_list field for bad sample mask")
    local l1b_hdf_file = self.config:l1b_hdf_file()
    local sid = self.config:l1b_sid_list()
    local sounding_num = sid:sounding_number()
    local bad_sample_mask = l1b_hdf_file:read_double_3d("/InstrumentHeader/bad_sample_list")(Range.all(), sounding_num, Range.all())
    
    return bad_sample_mask
end

------------------------------------------------------------
--- Retrieves a bad pixel mask out of an extra dimension
--- in the snr_coef dataset. This was how things were done
--- pre B8.00
------------------------------------------------------------

function OcoConfig.snr_coef_bad_sample_mask(self)
    local l1b_hdf_file = self.config:l1b_hdf_file()
    local sid = self.config:l1b_sid_list()
    local sounding_num = sid:sounding_number()
    local snr_coef = l1b_hdf_file:read_double_4d_sounding("/InstrumentHeader/snr_coef", sounding_num)

    if (snr_coef:depth() > 2) then
        self.config:diagnostic_message("Using snr_coeff third component for bad sample data")
        local bad_sample_mask = snr_coef(Range.all(), Range.all(), 2)
        return bad_sample_mask
    end
    
    self.config:diagnostic_message("No bad sample data found in snr_coeff")
    return nil
end

------------------------------------------------------------
--- Extends the snr coef bad sample mask routine but also
--- applies South Atlantic Anomoly (SAA) bad sample
--- removal.
---
--- This first applies the mask "bad_sample_mask_before_saa",
--- then it applies the SAA filtering, but only if we are
--- within the latitude/longitude box given by latitude_min, etc.
--- We use the given saa_tolerance to find bad samples.
------------------------------------------------------------

function OcoConfig.bad_sample_mask_with_saa(self)
   -- Ensure the saa_tolerance is set
   if (not self.saa_tolerance) then
      error("SAA filtering saa_tolerance not set")
   end

   -- Shape: band x sample_index
   local bad_samp_mask = self:bad_sample_mask_before_saa()

   if (bad_samp_mask ~= nil) then
      -- We haven't read in the l1b object yet, because we have
      -- a chicken an egg problem where we haven't yet generated
      -- the noise model. But go ahead a create a local copy that we 
      -- use just to get the latitude and longitude
      local sid = self.config:l1b_sid_list()
      local l1b_temp = Level1bOco(self.config:l1b_hdf_file(), sid)
      local lat = l1b_temp:latitude(0).value
      local lon = l1b_temp:longitude(0).value
      if(lat >= self.latitude_min and lat <= self.latitude_max and
	 lon >= self.longitude_min and lon <= self.longitude_max) then
	 self.config:diagnostic_message("Applying SAA bad sample filtering with saa_tolerance: " .. self.saa_tolerance)
	 
	 local frame_num = sid:frame_number()
	 local sounding_num = sid:sounding_number()
	 
	 for band_idx = 0, bad_samp_mask:rows() - 1 do
	    if l1b_temp:has_spike_eof(band_idx) then
	       local saa_weighted = l1b_temp:spike_eof(band_idx)
	       for samp_idx = 0, bad_samp_mask:cols() - 1 do
		  if(saa_weighted(samp_idx) >= self.saa_tolerance) then
		     bad_samp_mask:set(band_idx, samp_idx, 1.0)
		     self.config:diagnostic_message("Setting SAA bad point band_idx " .. band_idx .. " samp_idx " .. samp_idx)
		  end
	       end
	    else
	       self.config:diagnostic_message("Not applying SAA filtering for " .. band_idx .. "dataset does not exist")
	    end
	 end
      else
	 self.config:diagnostic_message("Not applying SAA bad sample filtering because we are outside of the lat/lon box we apply this to")
      end
   else
      self.config:diagnostic_message("Not applying SAA bad sample filtering due to nil bad sample list")
   end
   
   return bad_samp_mask
end

------------------------------------------------------------
--- Retrieve offset scaling for static HDF file
------------------------------------------------------------

function OcoConfig.oco_offset_scaling(field)
   return 
      function (self)
         local sid = self.config:l1b_sid_list()
         sounding_pos = sid:sounding_position()
         scaling_all = self.config:h():read_double_with_unit_2d(field)
         return ArrayWithUnit_1d(scaling_all.data(sounding_pos-1, Range.all()), scaling_all.units)
      end
end

------------------------------------------------------------
--- Determine if we are land or water. Return "land" or 
--- "water"
------------------------------------------------------------

function OcoConfig:land_or_water()
   local lf
   if (self.unscaled_l1b ~= nil) then
      lf = self.unscaled_l1b:land_fraction()
   else
      lf = self.l1b:land_fraction()
   end

   if(lf > 80.0) then
      return "land"
   elseif(lf < 20.0) then
      return "water"
   else
      error("Invalid land fraction value read from L1B file: " .. lf .. " (must be > 80 or < 20 to process Level 2)")
   end
end

------------------------------------------------------------
--- Defines the ground type name we should use, based on
--- land water fag
------------------------------------------------------------

function OcoConfig:ground_type_name()
   if(self:land_or_water() == "land") then
      return "brdf_soil"
   else
      return "coxmunk"
   end
end

------------------------------------------------------------
--- Fluorescence only for land runs
------------------------------------------------------------

OcoConfig.fluorescence_effect_land_only = ConfigCommon.fluorescence_effect:new()
function OcoConfig.fluorescence_effect_land_only:create()
   -- Only use this for land runs
   if(self.config:land_or_water() == "land") then
      return ConfigCommon.fluorescence_effect.create(self)
   else
      return nil
   end
end

------------------------------------------------------------
--- ZLO only for land runs
------------------------------------------------------------

OcoConfig.zero_offset_waveform_land_only = ConfigCommon.zero_offset_waveform:new()
function OcoConfig.zero_offset_waveform_land_only:create()
   -- Only use this for land runs
   if(self.config:land_or_water() == "land") then
      return ConfigCommon.zero_offset_waveform.create(self)
   else
      return nil
   end
end

------------------------------------------------------------
--- Create ground based on the surface type
------------------------------------------------------------

OcoConfig.ground_from_ground_type = DispatchCreator:new()

function OcoConfig.ground_from_ground_type:get_creator()
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

OcoConfig.ground_from_ground_type_scaled = DispatchCreator:new()

function OcoConfig.ground_from_ground_type_scaled:get_creator()
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
         return ConfigCommon.ground_coxmunk_scaled
      end
   else
      error("Invalid ground type value: " .. ground_type)
   end
end

------------------------------------------------------------
--- Create lambertian ground initial guess from radiance
--- but the other parts from the static HDF file
------------------------------------------------------------

function OcoConfig.oco_albedo_from_radiance(polynomial_degree)
   return function(self, band_idx)

      -- Fsun0 - Roughly at 1 AU
      local stokes_coef = self.config.l1b:stokes_coef()
      local solar_strength = {}
      for idx = 1, self.config.number_pixel:rows() do
         -- Create SolarDopplerShiftPolynomial so we can compute solar distance
         solar_doppler_shift = 
            SolarDopplerShiftPolynomial.create_from_l1b(self.config.l1b, idx-1, true)
         solar_dist = solar_doppler_shift:solar_distance().value

         -- Account for solar distance Fsun = Fsun0 / (solar_distance_meters/AU)^2
         solar_strength[idx] = self.config.fm.atmosphere.ground.solar_strength[idx]
         solar_strength[idx] = solar_strength[idx] / solar_dist^2
 
         -- Account for stokes element for I
         solar_strength[idx] = solar_strength[idx] * stokes_coef(idx-1, 0)
      end

      local continuum_points = {
         -- ABO2
         { { Range(20,31) }, { Range(912,923) } },
         -- WCO2
         { { 56, Range(85,88), Range(99,104), 145, Range(151,158) },
           { Range(766,770), 892, 893, Range(957,963), Range(967,969), Range(1000,1002) } },
         -- SCO2
         { { 10, 11, 12, 25, 26, 86, 87, 95, 105, 127, 128, 129 },
           { 949, 950, Range(973,976), 992, 993, 997, 998, 1008, 1009 } },
      }
      local use_range_max = { false, false, false }

      return ConfigCommon.albedo_from_radiance(self, band_idx, solar_strength, continuum_points, use_range_max, polynomial_degree)
   end
end

------------------------------------------------------------
--- The number of dispersion values for OCO-2 is 6, but
--- for OCO-1 there were 10. This number is carried into
--- many older files. For backwards compatibility we
--- extend the dispersion covariance for these older files.
------------------------------------------------------------

function OcoConfig.dispersion_covariance_i(field)
    return function(self, i)
        local num_disp = self.config.l1b:spectral_coefficient():cols()
        local orig_disp = ConfigCommon.hdf_covariance_i(field)(self, i)
        if num_disp == orig_disp:rows() then
            return orig_disp
        else
            local new_disp = Blitz_double_array_2d(num_disp, num_disp)
            new_disp:set(Range.all(), Range.all(), 0.0)
            local orig_vals_range = Range(0, orig_disp:rows()-1)
            new_disp:set(orig_vals_range, orig_vals_range, orig_disp(Range.all(), Range.all()))
            return new_disp
        end
    end
end

------------------------------------------------------------
--- Create ILS table by reading from an L1B hdf file using
--- Lua to do the reading and index massaging.
------------------------------------------------------------
OcoConfig.ils_table_l1b = Creator:new()

function OcoConfig.ils_table_l1b:create()
   -- Load reference to L1B HDF file
   local l1b_hdf_file = self.config:l1b_hdf_file()

   local sid = self.config:l1b_sid_list()
   sounding_num = sid:sounding_number()

   -- Use a simple table ILS since the ILS supplies
   -- a full set of delta_lambda/repsonses per pixel
   interpolate = false

   local idx
   local wavenumber, delta_lambda
   local res = {}
   local scale = Blitz_double_array_1d(1)
   local scale_flag = Blitz_bool_array_1d(1)
   for idx = 0, self.config.number_pixel:rows() - 1 do
      scale:set(0, self.scale_apriori[idx+1])
      scale_flag:set(0, self.retrieve_bands[idx+1])

      local hdf_band_name = self.config.common.hdf_band_name:value(idx)
      local desc_band_name = self.config.common.desc_band_name:value(idx)

      -- Delta lambda
      delta_lambda = l1b_hdf_file:read_double_4d_sounding("/InstrumentHeader/ils_delta_lambda", sounding_num)
      delta_lambda = delta_lambda(idx, Range.all(), Range.all())

      -- Load ILS response values
      response = l1b_hdf_file:read_double_4d_sounding("/InstrumentHeader/ils_relative_response", sounding_num)
      response = response(idx, Range.all(), Range.all())

      -- Calculate center wavenumber from dispersion, there should be number pixel
      -- of these per spectrometer
      wavenumber = self.config.dispersion[idx+1]:pixel_grid():data()

--    res[idx+1] = IlsTableLinear(wavenumber, delta_lambda, response, scale, scale_flag, desc_band_name, hdf_band_name, interpolate)
      res[idx+1] = IlsTableLog   (wavenumber, delta_lambda, response, scale, scale_flag, desc_band_name, hdf_band_name, interpolate)
   end
   return res
end

-- For IlsInstrument, we need to build up the initial guess in
-- the same order things are put into the statevector.
function OcoConfig.ils_table_l1b:initial_guess_i(i)
   local res = CompositeInitialGuess()
   local scale = Blitz_double_array_1d(1)
   scale:set(0, self.scale_apriori[i])
   local scale_flag = Blitz_bool_array_1d(1)
   scale_flag:set(0, self.retrieve_bands[i])
   local covariance = Blitz_double_array_2d(1,1)
   covariance:set(0, 0, self.scale_cov[i])
   local ig = InitialGuessValue()
   ig:apriori_subset(scale_flag, scale)
   ig:apriori_covariance_subset(scale_flag, covariance)
   res:add_builder(ig)
   return res
end

function OcoConfig.ils_table_l1b:register_output(ro)
   if(self.use_scale) then
      for i = 1, self.config.number_pixel:rows() do
	 local hdf_band_name = self.config.common.hdf_band_name:value(i-1)
	 ro:push_back(IlsTableLogOutput.create(self.config.ils_func[i], hdf_band_name))
      end
   end
end

------------------------------------------------------------
--- Open a ECWMF based on the orbit simulator meteorology 
--- file. 
------------------------------------------------------------

OcoConfig.oco_meteorology = Creator:new()

function OcoConfig.oco_meteorology:create()
   local sid = self.config:l1b_sid_list()
   if (self.config.met_file) then
       local met = OcoSimMetEcmwf(self.config.met_file, self.config:l1b_sid_list())
       self.config.input_file_description = self.config.input_file_description .. 
          "Meteorology input file:    " .. self.config.met_file .. "\n"
       return met
   end
end

function OcoConfig.oco_meteorology:register_output(ro)
    if (self.config.met) then
        ro:push_back(MetPassThroughOutput(self.config.met))
    end
end

------------------------------------------------------------
--- Get the CO2 apriori from the profile given in the scene
--- file, rather than trying to guess it. This is useful 
--- when we are trying to exactly match a OCO simulator run.
---
--- If you use this function, then you should have scene_file
--- defined in the config object
------------------------------------------------------------

function OcoConfig:co2_apriori_from_scene()
   local t = OcoSimApriori(self.config.scene_file, self.config:l1b_sid_list())
   return t:co2_vmr_grid(self.config.pressure)
end

------------------------------------------------------------
--- Use tropopause height for initial guess as cirrus ice
--- height
------------------------------------------------------------

function OcoConfig.tropopause_height_ap(self, base, type, aer_name)
    local ap = self:apriori_initial(base, type, aer_name)
    
    local t = self.config:reference_co2_apriori_met_obj()
    local psurf = self.config.met:surface_pressure()
    ap:set(1, (t:tropopause_pressure() - 100) / psurf)
    
    return ap
end


