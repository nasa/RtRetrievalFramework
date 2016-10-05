------------------------------------------------------------
--- Loads state apriori from a previous L2 retrieval output
--- file.
------------------------------------------------------------

require "gosat_base_config"

config = GosatBaseConfig:new()

-- Turn on diagnostics for this specicial type of config
config.diagnostic = true

-- Load retrieval results from L2 output file
l2_retrieval_file = os.getenv("l2_retrieval_file")

if not l2_retrieval_file then
   error('The "l2_retrieval_file" environmental variable must be defined')
end

l2_ret = HdfFile(l2_retrieval_file)
ret_sv_names  = l2_ret:read_string_vector("/RetrievedStateVector/state_vector_names", 0)
ret_sv_values = l2_ret:read_double_2d("/RetrievedStateVector/state_vector_result")

if(config.diagnostic) then
   print("State Vector Names:")
   print(ret_sv_names)
   print("Retrieved State Vector Values:")
   print(ret_sv_values)
end

-------------------
-- Generic routines
-------------------

function get_sv_values_by_name(search_name)
   -- Find matching sv indexes
   sv_indexes = {}
   for sv_idx = 0, ret_sv_names:size()-1 do
      curr_name = ret_sv_names:value(sv_idx)
      if curr_name == search_name or curr_name:find(search_name, 1, true) then
	 table.insert(sv_indexes, sv_idx)
      end
   end

   -- Build result array
   sv_sub = Blitz_double_array_1d(#sv_indexes)
   for sub_idx, sv_idx in ipairs(sv_indexes) do
      sv_sub:set(sub_idx-1, ret_sv_values(0, sv_idx))
   end

   return sv_sub
end

function retrieved_ap(sv_name, orig_ap_func, range)
   return function(self)
      local ap_vals
      if orig_ap_func then
	 ap_vals = orig_ap_func(self)
      end
      local ret_vals = get_sv_values_by_name(sv_name)
      
      if ret_vals:rows() > 0 then
	 if(ap_vals and self.config.diagnostic) then
	    print(sv_name .. " Original:")
	    print(ap_vals)
	 end
      
	 if ap_vals then
	    if range then
	       ap_vals:set(range, ret_vals)
	    else
	       ap_vals:set(Range.all(), ret_vals)
	    end
	 else
	    ap_vals = ret_vals
	 end

	 if(self.config.diagnostic) then
	    print(sv_name .. " Updated:")
	    print(ap_vals, "\n")
	 end
      end

      if not ap_vals then
	 error("Attempting to update " .. sv_name .. " with a nil value")
      end

      return ap_vals
   end
end

function retrieved_ap_i(sv_name, orig_ap_func, range)
   return function(self, i)
      local desc_band_name = self.config.common.desc_band_name:value(i)
      local orig_ap_func_i
      if orig_ap_func then
	 orig_ap_func_i = function(self)
	    return orig_ap_func(self, i)
	 end
      end
      return retrieved_ap(sv_name .. " " .. desc_band_name, orig_ap_func_i, range)(self)
   end
end

function retrieved_gsd(sv_name, orig_gsd_func, post_str)
   if not post_str then
      post_str = ""
   end
   return function(self)
      local gsd_obj = orig_gsd_func(self)
      print("gsd obj (" .. sv_name..")= ", gsd_obj)
      local gsd_ig = gsd_obj:initial_guess()
      local gsd_ap = gsd_obj:apriori()
      for j=0, self.config.number_pixel:rows()-1 do
	 local desc_band_name = self.config.common.desc_band_name:value(j)
	 local band_vals = get_sv_values_by_name(sv_name .. " " .. desc_band_name .. post_str)
	 if band_vals:rows() > 0 then
	    if(self.config.diagnostic) then
	       print(sv_name .. " " .. j .. " Original: ")
	       print(gsd_ig(j, 0, Range.all()))
	    end
	    
	    gsd_ig:set(j, 0, Range.all(), band_vals)
	    gsd_ap:set(j, 0, Range.all(), band_vals)
	    
	    if(self.config.diagnostic) then
	       print(sv_name .. " " .. j .. " Updated: ")
	       print(gsd_ig(j, 0, Range.all()), "\n")
	    end
	 end
      end
      
      gsd_obj:initial_guess(gsd_ig)
      gsd_obj:apriori(gsd_ap)
      
      return gsd_obj
   end
end

------------------------
-- Specific values setup
------------------------

-- Set Zero level offset to retrieved value
config.fm.instrument.instrument_correction.zero_offset_waveform.apriori = retrieved_ap_i("Zero offset waveform scale factor", config.fm.instrument.instrument_correction.zero_offset_waveform.apriori)

-- Set surface pressure to retrieved value
config.fm.atmosphere.pressure.apriori = retrieved_ap("Surface Pressure")

-- Set temperature offset to retrieved value
config.fm.atmosphere.temperature.apriori = retrieved_ap("Temperature Offset")

-- Set CO2 levels to retrieved values
config.fm.atmosphere.absorber.CO2.apriori = retrieved_ap("CO2 VMR")

-- Set H2O scaling to retrieved value
sv_name = "H2O Scaling factor"
config.fm.atmosphere.absorber.H2O.scale_apriori = get_sv_values_by_name(sv_name)(0)
if config.diagnostic then
   print(sv_name .. " Updated:")
   print(config.fm.atmosphere.absorber.H2O.scale_apriori, "\n")
end

-- Set aerosol shape values
for aer_idx, aer_name in ipairs(config.fm.atmosphere.aerosol.aerosols) do
   config.fm.atmosphere.aerosol[aer_name].apriori = retrieved_ap("Aerosol Shape " .. aer_name)
end

-- Set dispersion to retrieved value
orig_disp_ap_func = config.fm.instrument.dispersion.apriori
function retrieved_dispersion_ap(self, i)
   if(not self.retrieved_disp) then
      sv_name = "Instrument Dispersion"

      local disp_vals = orig_disp_ap_func(self)
      if(self.config.diagnostic) then
	 print(sv_name .. " Original:")
	 print(disp_vals.data)
      end

      for j=0, self.config.number_pixel:rows()-1 do
	 local desc_band_name = self.config.common.desc_band_name:value(j)
	 local band_disp_offset = get_sv_values_by_name(sv_name .. " " .. desc_band_name)
	 if band_disp_offset:rows() == 0 then
	    error("Could not find state vector item: " .. sv_name .. " " .. desc_band_name)
	 end
	 disp_vals.data:set(j, 0, band_disp_offset(0))
      end
      
      if(self.config.diagnostic) then
	 print(sv_name .. " Updated:")
	 print(disp_vals.data, "\n")
      end
      self.retrieved_disp = disp_vals
   end
   return ArrayWithUnit_1d(self.retrieved_disp.data(i, Range.all()), self.retrieved_disp.units)
end
config.fm.instrument.dispersion.apriori = retrieved_dispersion_ap
config.fm.instrument.dispersion.creator = ConfigCommon.dispersion_polynomial

-- Set ground spectrally dependent values 
config.fm.atmosphere.ground.lambertian_gsd = retrieved_gsd("Ground Lambertian", config.fm.atmosphere.ground.lambertian_gsd)
config.fm.atmosphere.ground.coxmunk_gsd = retrieved_gsd("Ground Cox-Munk", config.fm.atmosphere.ground.coxmunk_gsd, " Lambertian Albedo")
config.fm.atmosphere.ground.coxmunk_plus_lamb_gsd = retrieved_gsd("Ground Cox-Munk", config.fm.atmosphere.ground.coxmunk_plus_lamb_gsd, " Lambertian Albedo")

-- Set windspeed from retrieved value
config.fm.atmosphere.ground.windspeed.apriori = retrieved_ap("Ground Cox-Munk Windspeed", config.fm.atmosphere.ground.windspeed.apriori)

config:do_config()

print("Updated State Vector:", config.state_vector)