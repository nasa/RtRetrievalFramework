------------------------------------------------------------
--- Uses the orbit simulator base configuration instead
--- of the production base configuration.
------------------------------------------------------------

require "oco_baseline_config"

config = OcoBaselineConfig:new()

config.diagnostic = true
config.solver_constant.max_iteration = 1
table.insert(config.fm.instrument.instrument_correction.ic, "radiance_scaling")
config.fm.atmosphere.aerosol.creator = ConfigCommon.rayleigh_only

config.scene_file = os.getenv("scene_file")
config.albedo_file = os.getenv("albedo_file")

config.write_jacobian = false
config.do_retrieval = false

function config:scene()
   if(self.scene_file and not self.scene_v) then 
      self.scene_v = HdfFile(self.scene_file) 
   end
   return self.scene_v
end

function config:albedo()
   if(self.albedo_file and not self.albedo_v) then 
      self.albedo_v = HdfFile(self.albedo_file) 
   end
   return self.albedo_v
end

function config:os_exposure_index()
   if(not self.exposure_index) then
      local sid = self:l1b_sid_list()
      local l1b_hdf = self:l1b_hdf_file()
      num_snd = l1b_hdf:read_double_2d("/SoundingGeometry/sounding_id"):cols()
      self.exposure_index =  sid:frame_number()*num_snd + sid:sounding_number()
   end
   return self.exposure_index
end

function config:os_num_level()
   if(not self.num_level) then
      local scene = self:scene()
      self.num_level = scene:read_int_1d("Simulation/Thermodynamic/num_layers")(self:os_exposure_index())+1
   end
   return self.num_level
end

function config:os_scene_pressure()
   local scene = self.config:scene()
   local press_lev = scene:read_double_2d("Simulation/Thermodynamic/pressure_level")(self.config:os_exposure_index(), Range(1,self.config:os_num_level()-1))
   print("OS Pressure:", press_lev)
   return press_lev
end

function config:os_surface_pressure()
   nlev = self.config:os_num_level()-1
   psurf = config.os_scene_pressure(self)(Range(nlev-1,nlev-1))
   print("OS Surface pressure:", psurf)
   return psurf
end

function config:os_temperature()
   local scene = self.config:scene()
   local temp_lev = scene:read_double_2d("Simulation/Thermodynamic/temperature_level")(self.config:os_exposure_index(), Range(1,self.config:os_num_level()-1))
   print("OS Temperature:", temp_lev)
   return temp_lev
end

function config.os_gas(gas_name)
   return function(self)

       if(not self.os_gas_table) then
	  local scene = self.config:scene()
	  local gas_names = scene:read_string_vector("Simulation/Gas/species_id", self.config:os_exposure_index())
	  local gas_data = scene:read_double_3d("Simulation/Gas/species_density")
	  self.os_gas_table = {}
	  for gidx = 0, gas_names:size()-1 do
	     curr_name = gas_names:value(gidx):gsub(" ","")
	     self.os_gas_table[curr_name] = gas_data(self.config:os_exposure_index(), gidx, Range(0,self.config:os_num_level()-2))
	  end
    	  for curr_name, gas_vals in pairs(self.os_gas_table) do
	     if(curr_name:find("AIR") == nil) then
		gas_vals = gas_vals / self.os_gas_table["AIR_dry"]
		self.os_gas_table[curr_name] = gas_vals
	     end
	  end
       end
       print("OS Gas ", gas_name, self.os_gas_table[gas_name])
       return self.os_gas_table[gas_name]

   end
end   

function config:os_windspeed()
   local scene = self.config:scene()
   local exposure_idx = self.config:os_exposure_index()
   return scene:read_double_1d("/Simulation/Surface/wind_speed")(Range(exposure_idx, exposure_idx))
end


function config.os_lambertian(field)
   return function(self)
       local alb = self.config:albedo()
       local sid = self.config:l1b_sid_list()

       orig_gsd = OcoConfig.oco_ground_from_radiance(field)(self)
       print("Calculated albedo:", orig_gsd:apriori())

       local ap = alb:read_double_4d("/Simulation/Surface/albedo_polynomial")(sid:frame_number(), sid:sounding_number(), Range.all(), Range.all())
       local eff = alb:read_double_4d("/Simulation/Surface/eff_albedo")(sid:frame_number(), sid:sounding_number(), Range.all(), Range.all())
       local ig = ap
       print("Read albedo:", ap)
       print("Effective albedo:", eff)

       local band_ref = self.config.common.band_reference

       local apcov = self.config:h():read_double_3d(field .. "/covariance")

       return GroundSpectrallyDependentParameter(ap, ig, apcov, band_ref, 1, 2)
   end
end


config.fm.atmosphere.pressure = {
   pressure_levels = config.os_scene_pressure,
   apriori = config.os_surface_pressure,
   covariance = ConfigCommon.hdf_covariance("Surface_Pressure"),
   creator = ConfigCommon.pressure_fixed_level,
}

config.fm.atmosphere.temperature = {
   creator = ConfigCommon.temperature_fixed_level,
   levels = { apriori= config.os_temperature, },
   offset = { apriori =ConfigCommon.hdf_apriori("Temperature/Offset"),
	      covariance = ConfigCommon.hdf_covariance("Temperature/Offset"), },
}

config.fm.atmosphere.absorber = {
   creator = ConfigCommon.absorber_creator,
   gases = {"CO2", "H2O", "O2"},
   CO2 = {
      apriori = config.os_gas("CO2"),
      covariance = ConfigCommon.hdf_covariance("Gas/CO2"),
      absco = "v3.3.0/lowres/co2_v3.3.0-lowres.hdf",
      creator = ConfigCommon.vmr_level,
   },
   H2O = {
      apriori = config.os_gas("H2O"),
      scale_apriori = 1.0,
      scale_cov = 0.25,
      absco = "v3.3.0/lowres/h2o_v3.3.0-lowres.hdf",
      creator = ConfigCommon.vmr_fixed_level_scaled,
   },
   O2 = {      
      apriori = config.os_gas("O2"),
      absco = "v3.3.0/lowres/o2_v3.3.0-lowres.hdf",
      creator = ConfigCommon.vmr_fixed_level_constant_well_mixed,
   },
}

config.fm.atmosphere.ground.windspeed.apriori = config.os_windspeed
--config.fm.atmosphere.ground.lambertian_gsd = config.os_lambertian("Ground/Lambertian")

config:do_config()
