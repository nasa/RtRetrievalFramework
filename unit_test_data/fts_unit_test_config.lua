-- This adds FTS specific routines for use by base_config

require "config_common"

------------------------------------------------------------
--- Get the file names found in source_files, because
--- they change from run to run.
------------------------------------------------------------

FtsUnitTestConfig = ConfigCommon:new()

------------------------------------------------------------
--- Get spectral window from ASCII file
------------------------------------------------------------

FtsUnitTestConfig.fts_spectral_window = Creator:new()

function FtsUnitTestConfig.fts_spectral_window:create()
   local data = Blitz_double_array_3d(1,1,2)
   local win_ranges = ArrayWithUnit_3d(data, "cm^-1")
   win_ranges.value:set(0, 0, 0, 6166.0)
   win_ranges.value:set(0, 0, 1, 6286.0)
   return SpectralWindowRange(win_ranges)
end

------------------------------------------------------------
--- Create a level 1b file using the FTS spectra files
--- and runlog.
------------------------------------------------------------

FtsUnitTestConfig.level1b_fts = CreatorL1b:new()

function FtsUnitTestConfig.level1b_fts:sub_object_key()
   return {}
end

function FtsUnitTestConfig.level1b_fts:create_parent_object()
   local spectra_names = VectorString()
   for i,w in ipairs({self.config.spectrum_1_file, 
		      self.config.spectrum_2_file,
		      self.config.spectrum_3_file}) do
      spectra_names:push_back(w);
   end

   l1b_fts = Level1bFts(self.config.runlog_file, spectra_names, self.config.spec_win:spectral_bound());

   return l1b_fts
end

function FtsUnitTestConfig.level1b_fts:register_output(ro)
   ro:push_back(FtsRunLogOutput(self.config.l1b:run_log()))
end

------------------------------------------------------------
--- Number of pixels
------------------------------------------------------------

function FtsUnitTestConfig:l1b_number_pixel()
   local res = Blitz_int_array_1d(self.config.l1b:number_spectrometer())
   local i
   for i =1,self.config.l1b:number_spectrometer() do
      res:set(i - 1, self.config.l1b:radiance(i - 1):rows())
   end
   return res
end

------------------------------------------------------------
--- Create the FTS instrument.
------------------------------------------------------------

FtsUnitTestConfig.fts_instrument = CompositeCreator:new()

function FtsUnitTestConfig.fts_instrument:sub_object_key()
   return {"dispersion", "instrument_correction"}
end

function FtsUnitTestConfig.fts_instrument:create_parent_object(sub_object)
   local ils = VectorIls()
   for i =1,self.config.l1b:number_spectrometer() do
      ils:push_back(IlsFts.create(self.config.dispersion[i],
				  self.dispersion:perturb(),
				  self.config.l1b, i - 1,
				  self.config.common.desc_band_name:value(i-1),
				  self.config.common.hdf_band_name:value(i-1)))
   end
   return IlsInstrument(ils, self.config.instrument_correction)
end

function FtsUnitTestConfig.fts_instrument:add_to_statevector(sv)
   sv:add_observer(self.config.instrument)
   CompositeCreator.add_to_statevector(self, sv)
end

------------------------------------------------------------
--- Returns the surface pressure extracted from the runlog
------------------------------------------------------------

function FtsUnitTestConfig:surface_pressure_from_runlog()
   -- Convert mb to Pa
   local l1b = self.config.l1b
   local pout_arr

   -- Use a value only from a valid spectrometer
   for sidx = 0, l1b:number_spectrometer()-1 do
      if l1b:radiance(sidx):rows() > 0 then
	 pout = self.config.l1b:run_log(sidx).outside_pressure
	 pout_arr = Blitz_double_array_1d(1)
	 pout_arr:set(0, pout * 1e2)
	 break
      end
   end
   
   if pout_arr == nil then
      error("Could not find a valid spectometer for extracting surface pressure")
   end

   if (self.config.diagnostic) then
      print("Read surface pressure from FTS run log: ", pout_arr(0))
   end

   return pout_arr
end

FtsUnitTestConfig.chapman_boa_rt_sza_calculate = Creator:new()

function FtsUnitTestConfig.chapman_boa_rt_sza_calculate:create()
   if(self.config.diagnostic) then
      print("Solar zenith angles from L1B:")
      print(self.config.l1b:sza())
   end
   local sza = UplookingRaytracing.calculate_apparent_surface_sza(self.config.l1b, self.config.atmosphere)
   if(self.config.diagnostic) then
      print("Apparent surface solar zenith angles:")
      print(sza)
   end
   return ChapmanBoaRT.create(self.config.atmosphere, sza, self.config.spec_win:spectral_bound())
end

