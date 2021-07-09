--- There is nothing stopping you from just creating Lua object directly in
--- your config file. However there are pieces of code that get used in a
--- lot places, so we have some routines for common ways we build things.
--- You are free to use these, or just ignore them and do more complicated
--- things as you like.
---
--- This particular file contains a number of different routines, you need
--- to actually select which ones to use. See for example base_config.lua
--- which sets up the "standard" run, you can then create a config file
--- for any deviation from that.

------------------------------------------------------------
--- Class to hold stuff as we create it. When you create a instance of this
--- you should make sure to add anything used by one of the functions you
--- call below. Check each function for details.
------------------------------------------------------------

ConfigCommon = { diagnostic = false, input_file_description = "",
		 wrote_atmsphere_desc = false}

function ConfigCommon:new (o)
   o = o or {}   -- create object if user does not provide one
   setmetatable(o, self)
   self.__index = self
   return o
end

------------------------------------------------------------
-- Log diagnostic messages if requested.
------------------------------------------------------------

function ConfigCommon:diagnostic_message(str)
   if(self.diagnostic) then
      print("Lua diagnostic: " .. str)
      -- Can turn on memory usage, if we need to diagnose a problem
      -- print(io.open("/proc/self/status", "r"):read("*all"))
   end
end

------------------------------------------------------------
--- Utility function to find the directory for the current 
--- file.
------------------------------------------------------------

function ConfigCommon.local_dir(level)
  local level = level or 0
  local src = debug.getinfo(level + 2).source
  local res
  if(string.find(src, '/') == nil) then
     res = '.'
  else
     local src2 = string.gsub(src, '/[^/]+$', '')
     res = string.match(src2, '^@(.*)$')
  end
  return res
end

------------------------------------------------------------
--- We have lots of things that get created in "create a
--- object" and "do something to get the initial guess".
--- This captures that common behavior.
------------------------------------------------------------

Creator = {}

------------------------------------------------------------
--- Make a new creator. We should pass in the base ConfigCommon
--- class along with any table information needed by the
--- the object (e.g., "apriori")
------------------------------------------------------------

function Creator:new(o, config, name)
   o = o or {}   -- create object if user does not provide one
   setmetatable(o, self)
   self.__index = self
   o.config = config
   o.name = name
   return o
end

------------------------------------------------------------
--- Create the object
------------------------------------------------------------

function Creator:create()
   error("Need to define create in derived class")
end

------------------------------------------------------------
--- Return the initial guess. This defaults to a empty 
--- CompositeInitialGuess, which is the same thing as saying
--- we don't have an initial guess for this object (i.e., it
--- has no state vector elements in it)
------------------------------------------------------------

function Creator:initial_guess()
   return CompositeInitialGuess()
end

------------------------------------------------------------
-- Add anything we need to the give state vector.
------------------------------------------------------------

function Creator:add_to_statevector(sv)
end

------------------------------------------------------------
-- Register any output.
------------------------------------------------------------

function Creator:register_output(ro)
end

------------------------------------------------------------
--- It is a common pattern to have the particular choice
--- of a Creator depend on some dynamic property (e.g.,
--- land water fraction selecting ground type). This class
--- is a Creator which just uses a user supplied Creator to
--- dispatch each function.
------------------------------------------------------------
collectgarbage()
DispatchCreator = Creator:new()

function DispatchCreator:get_creator()
   error("Derived classes need to supply get_creator")
end

function DispatchCreator:table()
   return self:get_creator():new(self, self.config, self.name)
end

function DispatchCreator:create()
   return self:table():create()
end

function DispatchCreator:initial_guess()
   return self:table():initial_guess()
end

function DispatchCreator:add_to_statevector(sv)
   return self:table():add_to_statevector(sv)
end

function DispatchCreator:register_output(ro)
   return self:table():register_output(ro)
end

------------------------------------------------------------
--- It is a common pattern to have a object that is made up
--- of a number of other objects (e.g., Aerosol, made of 
--- AerosolExtinction and AerosolProperty). This class handles
--- that common case.
------------------------------------------------------------

CompositeCreator = Creator:new()

--- Function that lists keys to use for creating subobjects.
--- This is in the order it is created, if that matters for your
--- particular derived class.
function CompositeCreator:sub_object_key()
   error("Need to define sub_object_key in derived class")
end

--- Function that list keys in order we should put them in initial
--- guess. Default is just they same order as subobject_key, but this
--- gives a place to change the order if the build order and order in
--- the state vector aren't the same.
function CompositeCreator:sub_initial_guess_key()
   return self:sub_object_key()
end

--- Function to create the overall object once the sub objects
--- have been created.
function CompositeCreator:create_parent_object(sub_object)
   error("Need to define create_parent_object in derived class")
end

--- Create objects. We stash the sub objects in the self.config object using
--- the key, this allows things to be available for other classes to use
--- (e.g., create pressure and make available for other classes)
function CompositeCreator:create()
   local sub_object = {}, i, k
   for i,k in ipairs(self:sub_object_key()) do
      self.config:diagnostic_message("Creating " .. k)
      local t = self[k]
      if t ~= nil then
         local c = t.creator:new(t, self.config, k)
         sub_object[i] = c:create()
         self.config[k] = sub_object[i]
      end
   end
   return self:create_parent_object(sub_object)
end

function CompositeCreator:initial_guess()
   local res = CompositeInitialGuess()
   local i, k
   for i,k in ipairs(self:sub_initial_guess_key()) do
      self.config:diagnostic_message("Initial guess " .. k)
      local t = self[k]
      local c = t.creator:new(t, self.config, k)
      res:add_builder(c:initial_guess())
   end
   return res
end

function CompositeCreator:add_to_statevector(sv)
   local sub_object = {}, i, k
   for i,k in ipairs(self:sub_initial_guess_key()) do
      self.config:diagnostic_message("Adding to state vector " .. k)
      local t = self[k]
      local c = t.creator:new(t, self.config, k)
      c:add_to_statevector(sv)
   end
end

function CompositeCreator:register_output(ro)
   local sub_object = {}, i, k
   for i,k in ipairs(self:sub_object_key()) do
      self.config:diagnostic_message("Register output " .. k)
      local t = self[k]
      if t ~= nil then
         local c = t.creator:new(t, self.config, k)
         c:register_output(ro)
      end
   end
end

------------------------------------------------------------
--- This is a Creator where we get initial guess from a 
--- supplied aprior and covariance function. The default
--- is that we retrieve all the elements given by apriori.
------------------------------------------------------------

CreatorApriori = Creator:new { }

function CreatorApriori:retrieval_flag()
   local res
   if self.retrieved == nil or self.retrieved then
      res = ConfigCommon.create_flag(self:apriori_v():rows(), Range.all())
   else
      res = ConfigCommon.create_flag(self:apriori_v():rows())
   end
   return res
end

function CreatorApriori:iguess_v()
   if(self.iguess == nil) then
      return self:apriori_v()
   else
      return self:iguess()
   end
end

function CreatorApriori:apriori_v()
   return self:apriori()
end

function CreatorApriori:covariance_v()
   return self:covariance()
end

function CreatorApriori:initial_guess()
   local rflag = self:retrieval_flag()
   local ig = InitialGuessValue()
   ig:apriori_subset(rflag, self:apriori_v())
   ig:apriori_covariance_subset(rflag, self:covariance_v())
   ig:initial_guess_subset(rflag, self:iguess_v())
   return ig
end

------------------------------------------------------------
--- Creator for multiple spectrometer type objects
------------------------------------------------------------

CreatorMultiSpec = Creator:new()

function CreatorMultiSpec:retrieval_flag(i)
   flag = Blitz_bool_array_1d(self:apriori_v(i - 1):rows())

   -- Only retrieve this band if the configuration explicitly
   -- specifies that it should be retrieved. And if the
   -- particular band is non empty
   if self.retrieve_bands ~= nil and self.retrieve_bands[i] then
      flag:set(Range.all(), self.config:spec_flag()(i - 1))
   else
      flag:set(Range.all(), false)
   end

   return flag
end

function CreatorMultiSpec:iguess_v(i)
   if(self.iguess == nil) then
      return self:apriori_v(i)
   else
      return self:iguess(i)
   end
end

function CreatorMultiSpec:apriori_v(i)
   return self:apriori(i)
end

function CreatorMultiSpec:covariance_v(i)
   return self:covariance(i)
end

function CreatorMultiSpec:initial_guess(i)
   local ig
   if (i == nil) then
       -- Initial guess routine not called as part of a CompositeMultiSpec
       -- return all band's values
       ig = CompositeInitialGuess()
       for i=1,self.config.number_pixel:rows() do
           ig:add_builder(self:initial_guess(i))
       end
   else
      -- Called as part of a CompositeMultiSpec, return spectromteter's
      -- values
      local flag = self:retrieval_flag(i)
      ig = InitialGuessValue()
      ig:apriori_subset(flag, self:apriori_v(i - 1))
      ig:apriori_covariance_subset(flag, self:covariance_v(i - 1))
      ig:initial_guess_subset(flag, self:iguess_v(i - 1))
   end

   return ig 
end

------------------------------------------------------------
--- Overall function for doing configuration
------------------------------------------------------------

function ConfigCommon:do_config()
   self.register_output = VectorRegisterOutput()
   local t = self.fm
   local c = t.creator:new(t, self, "forward_model")

   self.forward_model = c:create()
   self.initial_guess = c:initial_guess()
   c:add_to_statevector(self.state_vector)
   self.state_vector:update_state(self.initial_guess:initial_guess(), 
                                  self.initial_guess:apriori_covariance())

   -- Setup forward model grid, based on initial state vector
   self.forward_model:setup_grid()

   if(self.do_retrieval) then
      self.solver:create(self)
      -- Only one kind of error_analysis right now, so just call directly.
      -- We can change this to pass in the configuration if needed.
      self:create_error_analysis()
   end

   c:register_output(self.register_output)

   -- To have the minimum requirements, L2FpConfigurationLua
   -- just requires a few global variables to be set. Copy
   -- our results to these variables.

   logger = FpLogger(self.log_level)
   self.forward_model:input_file_description(self.input_file_description)
   forward_model = self.forward_model
   solver = self.conn_solver
   -- This should be temporary, until we merge solver and iterative_solver
   iterative_solver = self.iter_solver
   stat_method_map = self.stat_method_map
   initial_guess = self.initial_guess
   number_pressure_level = self.number_pressure_level
   number_aerosol = self.number_aerosol
   number_band = self.spec_win:number_spectrometer()
   iteration_output = self.iteration_output
   register_output = self.register_output
end

------------------------------------------------------------
--- Open HDF file if we have one
------------------------------------------------------------

function ConfigCommon:h()
   if(self.static_file and not self.h_v) then 
      self.h_v = HdfFile(self.static_file)
      self.input_file_description = self.input_file_description .. 
	 "Static input file:   " .. self.static_file .. "\n"
   end   
   return self.h_v
end

------------------------------------------------------------
--- Open HDF file if we have one
------------------------------------------------------------

function ConfigCommon:h_solar()
   -- Use static_solar_file if found, otherwise use the same static input
   -- file that we use for everything else (self:h()).
   if(self.static_solar_file and not self.h_solar_v) then 
      self.h_solar_v = HdfFile(self.static_solar_file) 
      self.input_file_description = self.input_file_description .. 
	 "Solar input file:    " .. self.static_solar_file .. "\n"
   end
   if(not self.static_solar_file) then
      self.h_solar_v = self:h()
   end
   return self.h_solar_v
end

function ConfigCommon:h_eof()
   -- Use static_eof_file if found, otherwise use the same static input
   -- file that we use for everything else (self:h()).
   if(self.static_eof_file and not self.h_eof_v) then 
      self.h_eof_v = HdfFile(self.static_eof_file) 
      self.input_file_description = self.input_file_description .. 
	 "EOF input file:    " .. self.static_eof_file .. "\n"
   end
   if(not self.static_eof_file) then
      self.h_eof_v = self:h()
   end
   return self.h_eof_v
end

function ConfigCommon:h_imap()
   if(self.imap_file and self.imap_file ~= "" and not self.h_imap_v) then
      self.h_imap_v = HdfFile(self.imap_file)
      self.input_file_description = self.input_file_description .. 
        "IMAP input file:   " .. self.imap_file .. "\n"
   end
   return self.h_imap_v
end

function ConfigCommon:h_co2_profile()
   if(self.co2_pr_file and self.co2_pr_file ~= "" and not self.h_co2_profile_v) then
      self.h_co2_profile_v = HdfFile(self.co2_pr_file)
      self.input_file_description = self.input_file_description .. 
        "CO2 Profile input file:   " .. self.co2_pr_file .. "\n"
   end
   return self.h_co2_profile_v
end

function ConfigCommon:h_merra_aerosol()
   --- Use static_merra_aerosol_file if found, otherwise use h_aerosol()
   -- Use static_aerosol_file if found, otherwise use the same static input
   -- file that we use for everything else (self:h()).
   if(self.static_merra_aerosol_file and not self.h_merra_aerosol_v) then 
      self.h_merra_aerosol_v = HdfFile(self.static_merra_aerosol_file) 
      self.input_file_description = self.input_file_description .. 
	 "Merra Aerosol input file:  " .. self.static_merra_aerosol_file .. "\n"
   end
   if(not self.static_merra_aerosol_file) then
      self.h_merra_aerosol_v = self:h_aerosol()
   end
   return self.h_merra_aerosol_v
end

function ConfigCommon:h_aerosol()
   -- Use static_aerosol_file if found, otherwise use the same static input
   -- file that we use for everything else (self:h()).
   if(self.static_aerosol_file and not self.h_aerosol_v) then 
      self.h_aerosol_v = HdfFile(self.static_aerosol_file) 
      self.input_file_description = self.input_file_description .. 
	 "Aerosol input file:  " .. self.static_aerosol_file .. "\n"
   end
   if(not self.static_aerosol_file) then
      self.h_aerosol_v = self:h()
   end
   return self.h_aerosol_v
end

------------------------------------------------------------
-- Open Level 1b HDF file if we have one.
------------------------------------------------------------

function ConfigCommon:l1b_hdf_file()
   if(self.spectrum_file and not self.l1b_hdf_file_v) then
      self.l1b_hdf_file_v = HdfFile(self.spectrum_file)
      self.input_file_description = self.input_file_description .. 
	 "L1B input file:      " .. self.spectrum_file .. "\n"
   end
   return self.l1b_hdf_file_v
end

------------------------------------------------------------
-- Determine sounding id list from HDF file
-- Contents is instrument specific and not defined here
------------------------------------------------------------

function ConfigCommon:l1b_sid_list()
   error("l1b_sid_list definition is instrument specific and not defined in ConfigCommon.")
end

------------------------------------------------------------
-- Allows an apriori or covariance to be a function
-- that returns the value or a simple type
------------------------------------------------------------

function function_or_simple_value(val, ...)
    if type(val) == "function" then
        return val(...)
    else
        return val
    end
end

------------------------------------------------------------
--- Return various HDF creators
------------------------------------------------------------

function ConfigCommon.hdf_read_int_1d(field)
   return function(self)
             return self.config:h():read_int_1d(field)
          end
end

function ConfigCommon.hdf_read_double_scalar(field)
   return function(self)
             return self.config:h():read_double_scalar(field)
          end
end

function ConfigCommon.hdf_read_double_1d(field)
   return function(self)
             return self.config:h():read_double_1d(field)
          end
end

function ConfigCommon.hdf_read_spec_dom(field)
   return function(self)
             return SpectralDomain(self.config:h():read_double_with_unit_1d(field))
          end
end

function ConfigCommon.hdf_read_double_with_unit_1d(field)
   return function(self)
             return self.config:h():read_double_with_unit_1d(field)
          end
end

function ConfigCommon.hdf_read_double_with_unit(field)
   return function(self)
             return self.config:h():read_double_with_unit_1d(field)(0)
          end
end

function ConfigCommon.hdf_read_double_2d(field)
   return function(self)
             return self.config:h():read_double_2d(field)
          end
end

function ConfigCommon.hdf_read_double_3d(field)
   return function(self)
             return self.config:h():read_double_3d(field)
          end
end

function ConfigCommon.hdf_read_string(field)
   return function(self)
             return self.config:h():read_string(field)
          end
end

function ConfigCommon.hdf_read_string_vector(field)
   return function(self)
             local str_vec = self.config:h():read_string_vector(field)
             return str_vec
          end
end

function ConfigCommon.hdf_apriori_i(field)
   return function(self,i)
             return self.config:h():apriori(field, i) 
          end
end

function ConfigCommon.hdf_apriori_i_j(field, j)
   return function(self,i)
             return self.config:h():apriori(field, i, j) 
          end
end

function ConfigCommon.hdf_eof_apriori_i_j(field, j)
   return function(self,i)
             return self.config:h_eof():apriori(field, i, j) 
          end
end

function ConfigCommon.hdf_apriori(field)
   return function(self)
             return self.config:h():apriori(field) 
          end
end

function ConfigCommon.hdf_apriori_with_unit_i(field)
   return function(self,i)
             return self.config:h():apriori_with_unit(field, i) 
          end
end

function ConfigCommon.hdf_apriori_with_unit(field)
   return function(self)
             return self.config:h():apriori_with_unit(field) 
          end
end

function ConfigCommon.hdf_covariance_i(field)
   return function(self,i)
             return self.config:h():covariance(field, i) 
          end
end

function ConfigCommon.hdf_covariance_i_j(field, j)
   return function(self,i)
             return self.config:h():covariance(field, i, j) 
          end
end

function ConfigCommon.hdf_eof_covariance_i_j(field, j)
   return function(self,i)
             return self.config:h_eof():covariance(field, i, j) 
          end
end

function ConfigCommon.hdf_covariance(field)
   return function(self)
             return self.config:h():covariance(field) 
          end
end

function ConfigCommon.hdf_aerosol_apriori(base, type, aer_name)
   return function(self, config_name)
             if not aer_name then
                aer_name = config_name
             end

             if type then
                return self.config:h():apriori(base .. "/" .. aer_name .. "/" .. type)
             else
                return self.config:h():apriori(base .. "/" .. aer_name)
             end
          end
end

function ConfigCommon.ascii_aerosol_apriori(file)
   return function(self, a)
             return HeritageFile(file):data(a)
          end
end

function ConfigCommon.hdf_aerosol_covariance(base, type, aer_name)
   return function(self, config_name)
             if not aer_name then
                aer_name = config_name
             end

             if type then
                return self.config:h():covariance(base .. "/" .. aer_name .. "/" .. type)
             else
                return self.config:h():covariance(base .. "/" .. aer_name)
             end
          end
end

function ConfigCommon.hdf_aerosol_property(aer_group)
   return function(self, config_name)
             if not aer_group then
                aer_group = config_name
             end

             return AerosolPropertyHdf(self.config:h_aerosol(), aer_group .. "/Properties", self.config.pressure)
          end
end

------------------------------------------------------------
--- Variation of hdf_aerosol_property that allows the file
--- name to be passed in. Useful to one-shot sorts of tests
--- that we don't want to merge into our common aerosol file.
------------------------------------------------------------
function hdf_aerosol_property_file(fname, aer_group)
   return function(self, config_name)
             if not aer_group then
                aer_group = config_name
             end
	     local h = HdfFile(fname)
	     self.config.input_file_description = self.config.input_file_description ..
		"Aerosol input file:  " .. fname .. "\n"
             return AerosolPropertyHdf(h, aer_group .. "/Properties", self.config.pressure)
          end
end      

function ConfigCommon.ascii_dispersion(file, units)
   return function(self, i)
             data = HeritageFile(file):data("DISP_" .. (i + 1))
             return ArrayWithUnit_1d(data, units)
          end
end

function ConfigCommon.ascii_ground_apriori(file)
   return function(self, i)
             data = HeritageFile(file):data("ALBEDO_" .. (i + 1))
             return data
          end
end

function ConfigCommon.ascii_read(file, label, alt_label)
   return function(self)
             local ascii_file = HeritageFile(file)
         if ascii_file:column_index(label) >= 0 then
            return ascii_file:data(label)
         elseif alt_label ~= nil and ascii_file:column_index(alt_label) >= 0 then
            return ascii_file:data(alt_label)
         else
             if alt_label ~= nil then
                error("Could not find column named: " .. label .. " or alternative name: " .. alt_label .. " in file: " .. file)
             else 
                error("Could not find column named: " .. label .. " in file: " .. file)
             end
         end
          end
end

------------------------------------------------------------
-- Read from the atmosphere file defined in the configuration
------------------------------------------------------------

function ConfigCommon.read_atmosphere_file(field, alt_field_name)
   return function(self)
      if(self.config.diagnostic) then
         print( "Reading atmosphere file", self.config.atmosphere_file, field)
      end
      res = ConfigCommon.ascii_read(self.config.atmosphere_file, field, alt_field_name)(self)
      if(not self.config.wrote_atmsphere_desc) then
	 self.config.input_file_description = self.config.input_file_description .. 
	    "Atmosphere input file: " .. self.config.atmosphere_file .. "\n"
	 self.config.wrote_atmsphere_desc = true
      end
      return res
   end
end

------------------------------------------------------------
-- Read surface pressure from last level of atmophere
-- Pressure column
------------------------------------------------------------

function ConfigCommon.surface_pressure_from_atmosphere_file(field)
   return function(self)
      if not field then
         field = "Pressure"
      end

      if(self.config.diagnostic) then
         print( "Reading atmosphere file", self.config.atmosphere_file, field)
      end
      local atm_file = HeritageFile(self.config.atmosphere_file)
      if atm_file:has_value("Surface_Pressure") then
         local psurf = atm_file:value_double("Surface_Pressure")
         local psurf_arr = Blitz_double_array_1d(1)
         psurf_arr:set(0, psurf)
         return psurf_arr
      else
         local pressure = atm_file:data(field)

         -- Return as array of size 1
         psurf_idx = pressure:rows()-1
         return pressure(Range(psurf_idx,psurf_idx))
      end
   end
end

------------------------------------------------------------
--- Use CO2 from atmosphere file, or if not present, use
--- TCCON apriori
------------------------------------------------------------

function ConfigCommon:co2_from_atmosphere_or_tccon(co2_column_name)
    if (co2_column_name == nil) then
        co2_column_name = "CO2"
    end

    return function(self)
        ascii_obj = HeritageFile(self.config.atmosphere_file)
        if (ascii_obj:column_index(co2_column_name) >= 0) then
            return ConfigCommon.ascii_read(self.config.atmosphere_file, co2_column_name)(self)
        else
            return ConfigCommon.tccon_co2_apriori(self)
        end
    end
end

------------------------------------------------------------
-- Return an empty array
------------------------------------------------------------

function ConfigCommon.empty_array()
   return function(self)
             return Blitz_double_array_1d(0)
          end
end

------------------------------------------------------------
--- Create a blitz array from a lua array.
------------------------------------------------------------

function ConfigCommon.lua_to_blitz_double_1d(val)
   local nrow = #val
   local res = Blitz_double_array_1d(nrow)
   local i
   for i=1,nrow do
      res:set(i-1,val[i])
   end
   return res
end

------------------------------------------------------------
--- Create a blitz array from a lua array.
------------------------------------------------------------

function ConfigCommon.lua_to_blitz_double_2d(val)
   local nrow = #val
   local ncol = #val[1]
   local res = Blitz_double_array_2d(nrow, ncol)
   local i, j
   for i=1,nrow do
      for j=1,ncol do
	 res:set(i-1,j-1,val[i][j])
      end
   end
   return res
end

------------------------------------------------------------
--- Create a new InitialGuess that has the same apriori and
--- covariance as another one, but updates the first guess
------------------------------------------------------------

function ConfigCommon.update_initial_guess(ig, new_first_guess_value)
   local igb = InitialGuessValue()
   igb.apriori = ig:apriori()
   igb.initial_guess = new_first_guess_value
   igb.apriori_covariance = ig:apriori_covariance()
   local ignew = CompositeInitialGuess()
   ignew:add_builder(igb)
   return ignew
end

------------------------------------------------------------
--- Get various a priori values from ECMWF
------------------------------------------------------------

function ConfigCommon:met_pressure()
   local r = Blitz_double_array_1d(1)
   r:set(0,self.config.met:surface_pressure())
   return r
end

function ConfigCommon:met_windspeed()
   local r = Blitz_double_array_1d(1)
   r:set(0,self.config.met:windspeed())
   return r
end

function ConfigCommon:met_temperature()
   return self.config.met:temperature(self.config.pressure:pressure_level())
end

function ConfigCommon:met_h2o_vmr()
   return self.config.met:vmr("H2O", self.config.pressure:pressure_level())
end

-- These use the deprecated fixed pressure levels object
function ConfigCommon:met_temperature_fixed_pressure()
   return self.config.met:temperature(self.config.pinp:pressure_level())
end

function ConfigCommon:met_h2o_vmr_fixed_pressure()
   return self.config.met:vmr("H2O", self.config.pinp:pressure_level())
end

------------------------------------------------------------
--- Load the TCCON apriori object using ECMWF data
------------------------------------------------------------

function ConfigCommon:tccon_apriori_met()
    if(not self.tccon_ap_met_obj) then
        self.tccon_ap_met_obj = TcconApriori(self.met, self.l1b)
    end
    return self.tccon_ap_met_obj
end

------------------------------------------------------------
--- Load the TCCON apriori object using pressure/temp
--- objects.
------------------------------------------------------------

function ConfigCommon:tccon_apriori_pressure()
    if(not self.tccon_ap_met_obj) then
        self.tccon_ap_met_obj = TcconApriori(self.l1b, self.pressure, self.temperature)
    end
    return self.tccon_ap_met_obj
end

------------------------------------------------------------
--- Get tccon co2 apriori from tccon using ECMWF file
------------------------------------------------------------

function ConfigCommon:tccon_co2_apriori_met()
   local t = self.config:tccon_apriori_met()
   return t:co2_vmr_grid(self.config.pressure)
end

------------------------------------------------------------
--- Get tccon co2 apriori from tccon using pressure/temp objects
------------------------------------------------------------

function ConfigCommon:tccon_co2_apriori()
   local t = self.config:tccon_apriori_pressure()
   return t:co2_vmr_grid(self.config.pressure)
end

------------------------------------------------------------
--- Load the CO2 VMR object
------------------------------------------------------------

function ConfigCommon:reference_co2_apriori_met_obj()
    if (not self.ref_co2_ap_obj) then
        self.ref_co2_ap_obj = GasVmrApriori(self.met, self.l1b, self.altitude:value(0), self:h(), "/Reference_Atmosphere", "CO2")
    end

    return self.ref_co2_ap_obj 
end

------------------------------------------------------------
--- Get co2 apriori using reference apriori method
------------------------------------------------------------

function ConfigCommon:reference_co2_apriori_met_apriori()
   local t = self.config:reference_co2_apriori_met_obj()
   return t:apriori_vmr(self.config.pressure)
end

------------------------------------------------------------
--- Get co2 apriori for the profile file
------------------------------------------------------------

function ConfigCommon:co2_profile_file_apriori()
   local t = CO2ProfilePrior(self.config.met, self.config:h_co2_profile())
   return t:apriori_vmr(self.config.pressure)
end

------------------------------------------------------------
-- Short cut for making a flag array. This takes a size, creates
-- a bool array of that size and sets everything in the passed
-- in range to true, and everything else to false
------------------------------------------------------------

function ConfigCommon.create_flag(size, range_true)
   local res = Blitz_bool_array_1d(size)
   res:set(Range.all(), false)
   if range_true then
      res:set(range_true, true)
   end
   return res
end

------------------------------------------------------------
-- Flag for which spectrometers are valid with a non-empty
-- range
------------------------------------------------------------

function ConfigCommon:spec_flag()
   if not self.s_flag then
       if not self.spec_win then
          error("SpectralWindow not yet initialized, can not determing spectrometer flags")
       end

       local res = Blitz_bool_array_1d(self.spec_win:number_spectrometer())
       local sb = self.spec_win:spectral_bound()
       for spec_idx = 0, res:rows() -1 do
           range_size = sb:upper_bound(spec_idx).value - sb:lower_bound(spec_idx).value
           res:set(Range(spec_idx), range_size > 0)
       end
       self.s_flag = res
   end
   return self.s_flag
end

------------------------------------------------------------
--- Returns the value of the L1B object's spectral_coefficient
--- method.
------------------------------------------------------------

function ConfigCommon:l1b_spectral_coefficient()
   return self.config.l1b:spectral_coefficient_with_unit()
end

------------------------------------------------------------
--- Returns the value of the L1B object's spectral_coefficient
--- for a given spectrometer index.
------------------------------------------------------------

function ConfigCommon:l1b_spectral_coefficient_i(spec_idx)
   return self.config.l1b:spectral_coefficient_with_unit(spec_idx)
end

------------------------------------------------------------
--- Evaluates the table elements that are functions in the
--- context of the configuration environment.
--- This is useful for pull information out of the HDF
--- file and put it in a common location.
------------------------------------------------------------

ConfigCommon.table_function_eval = Creator:new()

function ConfigCommon.table_function_eval:create()
   result = {}
   for name, val in pairs(self) do
      if type(val) == "function" then
         result[name] = val(self)
      end
   end
   return result
end

------------------------------------------------------------
--- Create a constants class that is just the default constants 
------------------------------------------------------------

ConfigCommon.default_constant = Creator:new()

function ConfigCommon.default_constant:create()
   return DefaultConstant()
end

------------------------------------------------------------
--- Create ILS table by reading from hdf file.
------------------------------------------------------------
ConfigCommon.ils_table = Creator:new()

function ConfigCommon.ils_table:create()
   local i
   local res = {}
   local hdf_group = "Instrument/ILS"
   if(self.hdf_group_name) then
      hdf_group = self.hdf_group_name
   end
   for i=1,self.config.number_pixel:rows() do
      local desc_band_name = self.config.common.desc_band_name:value(i-1)
      local hdf_band_name = self.config.common.hdf_band_name:value(i-1)
      res[i] = IlsTableLinear(self.config:h(), i - 1, desc_band_name, hdf_band_name, hdf_group)
   end
   return res
end

function ConfigCommon.ils_table:initial_guess_i(i)
   return CompositeInitialGuess()
end

------------------------------------------------------------
--- Create dispersion as a polynomial
------------------------------------------------------------

ConfigCommon.dispersion_polynomial = Creator:new()

function ConfigCommon.dispersion_polynomial:retrieval_flag(i)
   -- Number of dispersion parameters
   local num_params = 1
   if self.num_parameters ~= nil then
       num_params = self.num_parameters
   end

   -- Retrieve only the first value for the dispersion, for all bands
   local ri_disp
   if self.retrieved == nil or self.retrieved then
      ri_disp = Range(0, num_params-1)
   else
      ri_disp = nil
   end

   -- Only retrieve dispersion when spectrometer range is non-empty
   local disp_coeff = self:coefficients(i)
   local disp_flag = ConfigCommon.create_flag(disp_coeff:rows(), ri_disp)
   if not self.config:spec_flag()(i - 1) then
      disp_flag:set(Range.all(), false)
   end  

   return disp_flag
end

function ConfigCommon.dispersion_polynomial:create()
   self.config.number_pixel = self:number_pixel()

   local i
   local res = {}
   for i=1,self.config.number_pixel:rows() do
      local disp_coeff = self:coefficients(i)
      local disp_units = self:units(i)
      local disp_flag = self:retrieval_flag(i)
      local desc_band_name = self.config.common.desc_band_name:value(i-1)
      res[i] = DispersionPolynomial(disp_coeff, disp_flag, disp_units,
                                    desc_band_name, 
                                    self.config.number_pixel(i - 1),
                                    self.is_one_based)
   end
   return res
end

function ConfigCommon.dispersion_polynomial:initial_guess()
   local i, j
   local res = CompositeInitialGuess()
   for i=1,self.config.number_pixel:rows() do
      local disp_coeff = self:coefficients(i)
      local disp_flag = self:retrieval_flag(i)
      local ig = InitialGuessValue()
      ig:apriori_subset(disp_flag, disp_coeff)
      ig:apriori_covariance_subset(disp_flag, self:covariance(i - 1))
      res:add_builder(ig)
   end
   return res
end

-- For IlsInstrument, we need to build up the initial guess in
-- the same order things are put into the statevector.
function ConfigCommon.dispersion_polynomial:initial_guess_i(i)
   local res = CompositeInitialGuess()
   local disp_coeff = self:coefficients(i)
   local disp_flag = self:retrieval_flag(i)
   local ig = InitialGuessValue()
   ig:apriori_subset(disp_flag, disp_coeff)
   ig:apriori_covariance_subset(disp_flag, self:covariance(i - 1))
   res:add_builder(ig)
   return res
end

function ConfigCommon.dispersion_polynomial:register_output(ro)
   for i=1,self.config.number_pixel:rows() do
      local hdf_band_name = self.config.common.hdf_band_name:value(i-1)
      ro:push_back(DispersionPolynomialOutput.create(self.config.dispersion[i], 
                                                     hdf_band_name))
   end
end

function ConfigCommon.dispersion_polynomial:coefficients(i)
   return self:apriori(i - 1).value
end

function ConfigCommon.dispersion_polynomial:units(i)
   return self:apriori(i - 1).units
end

------------------------------------------------------------
--- Create dispersion as a polynomial with instrument 
--- doppler applied to the offset term from the L1B
--- relative velocity.
------------------------------------------------------------

ConfigCommon.dispersion_polynomial_inst_doppler = ConfigCommon.dispersion_polynomial:new()

function ConfigCommon.dispersion_polynomial_inst_doppler:coefficients(i)
   local rel_vel = self.config.l1b:relative_velocity(0)
   local coeff = ConfigCommon.dispersion_polynomial.coefficients(self, i)
   local units = self:units(i)

   if(Unit(units):is_commensurate(Unit("cm^-1"))) then
      coeff:set(0, coeff(0) * (1 - rel_vel/Constants.speed_of_light))
   else
      coeff:set(0, coeff(0) / (1 - rel_vel/Constants.speed_of_light))
   end
   return coeff
end

------------------------------------------------------------
--- Create dispersion as a polynomial but fit the apriori
--- values according to solar lines in spectrum.
------------------------------------------------------------

ConfigCommon.dispersion_polynomial_fitted = ConfigCommon.dispersion_polynomial:new()

function ConfigCommon.dispersion_polynomial_fitted:coefficients(i)
   if(self.fitted_apriori == nil) then
      self.disp_fit = DispersionFit(self.config.l1b)
      ap_array = self:apriori().value
      self.fitted_apriori = self.disp_fit:fit(ap_array, 
                                              self.aband_solar_line_location, 
                                              self.aband_solar_line_width, 
                                              self.aband_search_width,
                                              self.aband_ils_offset,
                                              self.offset_scaling(self))
      if(self.config.diagnostic) then
         print("Fitting dispersion offset")
         print("------------------------")
         print("Original coefficients:")
         print(ap_array)
         print("Fitted coefficients:")
         print(self.fitted_apriori)
         print("Offset shifts:")
         print(self.disp_fit:shift())
         print("------------------------")
      end
   end
   return self.fitted_apriori(i - 1, Range.all())
end

function ConfigCommon.dispersion_polynomial_fitted:register_output(ro)
   ConfigCommon.dispersion_polynomial.register_output(self, ro)
   ro:push_back(DispersionFitOutput(self.disp_fit))
end

------------------------------------------------------------
--- Composites multiple CreatorMultiSpec created sub objects
------------------------------------------------------------

CompositeMultiSpecCreator = CompositeCreator:new()

function CompositeMultiSpecCreator:initial_guess()
   local res = CompositeInitialGuess()
   local i, k, s
   for s=1,self.config.l1b:number_spectrometer() do
       for i,k in ipairs(self:sub_initial_guess_key()) do
          self.config:diagnostic_message("Initial guess " .. k)
          local t = self[k]
          local c = t.creator:new(t, self.config, k)
          -- Initial guess may not be defined for all spectrometers
          local ig = c:initial_guess(s)
          if (ig) then
             res:add_builder(ig)
          end
       end
    end
   return res
end

----------------------------------------------------------
--- Create all the spectrum effects
------------------------------------------------------------

ConfigCommon.spectrum_effect_list = CompositeMultiSpecCreator:new()

function ConfigCommon.spectrum_effect_list:sub_object_key()
   return self.speceff
end

function ConfigCommon.spectrum_effect_list:create_parent_object(sub_object)

   local res = VectorVectorSpectrumEffect()
   local i
   for i=1,self.config.l1b:number_spectrometer() do
      local seff_vec = VectorSpectrumEffect()
      local j, s
      for j, s in ipairs(sub_object) do
         if s[i] then
            seff_vec:push_back(s[i])
         end
      end
      res:push_back(seff_vec)
   end
   return res
end

function ConfigCommon.spectrum_effect_list:add_to_statevector(sv)
   CompositeCreator.add_to_statevector(self, sv)
   sv:add_observer(self.config.spectrum_effect)
end

------------------------------------------------------------
--- Create all the instrument corrections.
------------------------------------------------------------

ConfigCommon.instrument_correction_list = CompositeMultiSpecCreator:new()

function ConfigCommon.instrument_correction_list:sub_object_key()
   return self.ic
end

function ConfigCommon.instrument_correction_list:create_parent_object(sub_object)
   
   local res = VectorVectorInstrumentCorrection()
   local i
   for i=1,self.config.number_pixel:rows() do
      local icorr_vec = VectorInstrumentCorrection()
      local j, s
      for j, s in ipairs(sub_object) do
	 if(s[i] ~= nil) then
	    icorr_vec:push_back(s[i])
	 end
      end
      res:push_back(icorr_vec)
   end
   return res
end

------------------------------------------------------------
--- Create the ILS instrument.
------------------------------------------------------------

ConfigCommon.ils_instrument = CompositeCreator:new()

function ConfigCommon.ils_instrument:sub_object_key()
   return {"dispersion", "ils_func", "instrument_correction"}
end

-- We handle dispersion and ils_func separately
function ConfigCommon.ils_instrument:sub_initial_guess_key()
   return {"instrument_correction"}
end

function ConfigCommon.ils_instrument:create_parent_object(sub_object)
   local ils = VectorIls()
   local i, ilf
   for i, ilf in ipairs(self.config.ils_func) do
      ils:push_back(IlsConvolution(self.config.dispersion[i], ilf, 
                                   self.ils_half_width[i]))
   end
   return IlsInstrument(ils, self.config.instrument_correction)
end

function ConfigCommon.ils_instrument:initial_guess()
   local res = CompositeInitialGuess()
   local i, j, k
   local t, c
   -- Reorder dispersion and ils_func. Because it comes in a IlsConvolution
   -- indexed by spectrometer
   for i=1,self.config.number_pixel:rows() do
      for j, k in ipairs({"dispersion", "ils_func"}) do
	 self.config:diagnostic_message("Initial guess " .. k .. "(" .. i .. ")")
	 t = self[k]
	 local c = t.creator:new(t, self.config, k)
	 res:add_builder(c:initial_guess_i(i))
      end
   end
   res:add_builder(CompositeCreator.initial_guess(self))
   return res
end

function ConfigCommon.ils_instrument:add_to_statevector(sv)
   CompositeCreator.add_to_statevector(self, sv)
   sv:add_observer(self.config.instrument)
end

------------------------------------------------------------
--- Create a ZeroOffsetWaveform correction
------------------------------------------------------------

ConfigCommon.zero_offset_waveform = CreatorMultiSpec:new()

function ConfigCommon.zero_offset_waveform:create()
   local res = {}
   for i=1,self.config.number_pixel:rows() do
      local fit_scale = self:retrieval_flag(i)
      res[i] = ZeroOffsetWaveform(self:apriori(i - 1)(0), fit_scale(0),
                                  self.config.dispersion[i], self.config:h(),
                                  i - 1, 
                                  self.config.common.desc_band_name:value(i-1),
                                  "Instrument/ZeroLevelOffset")
   end
   self.zero_offset_waveform = res
   return res
end

function ConfigCommon.zero_offset_waveform:initial_guess(i)
   if (not self.zero_offset_waveform) then
       return nil
   end
   return CreatorMultiSpec.initial_guess(self, i)
end

function ConfigCommon.zero_offset_waveform:register_output(ro)
   if(self.zero_offset_waveform) then
      for i=1,self.config.number_pixel:rows() do
	 ro:push_back(ZeroOffsetWaveformOutput.create(self.zero_offset_waveform[i], 
                           self.config.common.hdf_band_name:value(i-1)))
      end
   end
end

------------------------------------------------------------
--- Create a EmpiricalOrthogonalFunction correction
------------------------------------------------------------

ConfigCommon.empirical_orthogonal_function = CreatorMultiSpec:new()

function ConfigCommon.empirical_orthogonal_function:create()
   local res = {}
   local hdf_group
   local hdf_file

   if(self.hdf_group == nil) then
      hdf_group = "Instrument/EmpiricalOrthogonalFunction"
   else
      hdf_group = self.hdf_group
   end

   if(self.hdf_file == nil) then
      hdf_file = self.config:h_eof()
   else
      if type(self.hdf_file) == "function" then
         hdf_file = self.hdf_file(self)
      else
         hdf_file = self.hdf_file
      end
   end

   for i=1,self.config.number_pixel:rows() do
      local fit_scale = self:retrieval_flag(i)
      if(self.eof_used ~= nil and self.eof_used[i] == false) then
	 res[i] = nil
      elseif(self.by_pixel) then
	 if(self.scale_uncertainty) then
    -- Luabind can only handle up to 10 arguments per function. As an easy
    -- work around we put various values into an array
	    local mq = Blitz_double_array_1d(5)
	    mq:set(0, self:apriori(i - 1)(0))
	    mq:set(1, i - 1)
	    mq:set(2, self.config.sid:sounding_number())
	    mq:set(3, self.order)
	    mq:set(4, self.scale_to_stddev)
	    res[i] = EmpiricalOrthogonalFunction.create(
                                  fit_scale(0),
                                  hdf_file,
				  self.config.l1b:uncertainty_with_unit(i-1),
				  self.config.common.desc_band_name:value(i-1),
				  hdf_group,
				  mq)
	 else
	    res[i] = EmpiricalOrthogonalFunction(self:apriori(i - 1)(0),
                                  fit_scale(0),
                                  hdf_file,
                                  i - 1,
                                  self.config.sid:sounding_number(),
                                  self.order, 
                                  self.config.common.desc_band_name:value(i-1),
                                  hdf_group)
	 end
      else
         res[i] = EmpiricalOrthogonalFunction(self:apriori(i - 1)(0),
                                  fit_scale(0),
                                  self.config.dispersion[i],
                                  hdf_file,
                                  i - 1,
                                  self.config.sid:sounding_number(),
                                  self.order, 
                                  self.config.common.desc_band_name:value(i-1),
                                  hdf_group)
      end
   end
   self.eof = res
   return res
end

function ConfigCommon.empirical_orthogonal_function:register_output(ro)
   for i=1,self.config.number_pixel:rows() do
      if(self.eof[i] ~= nil) then
	 ro:push_back(EmpiricalOrthogonalFunctionOutput.create(self.eof[i], 
	              self.config.common.hdf_band_name:value(i-1)))
      end
   end
end

function ConfigCommon.empirical_orthogonal_function:retrieval_flag(i)
   --- Use CreatorMultiSpec
   local flag
   flag = CreatorMultiSpec.retrieval_flag(self, i)
   --- But override to make sure we don't retrieve EOFs we aren't using
   if(self.eof_used ~= nil and self.eof_used[i] == false) then
      flag:set(Range.all(), false)
   end
   return flag
end

------------------------------------------------------------
--- Create a RadianceScalingSvFit correction
------------------------------------------------------------

ConfigCommon.radiance_scaling_sv_fit = CreatorMultiSpec:new()

function ConfigCommon.radiance_scaling_sv_fit:create()
   local res = {}
   for i=1,self.config.number_pixel:rows() do
      local fit_scaling = self:retrieval_flag(i)
      local band_ref = self.config.common.band_reference
      res[i] = RadianceScalingSvFit(self:apriori(i - 1), 
                                     fit_scaling,
                                     band_ref(i-1),
                                     self.config.common.desc_band_name:value(i-1))
   end
   self.radiance_scaling = res
   return res
end

function ConfigCommon.radiance_scaling_sv_fit:register_output(ro)
    if (self.radiance_scaling) then
       for i=1,self.config.number_pixel:rows() do
          ro:push_back(RadianceScalingOutput.create(self.radiance_scaling[i], 
                       self.config.common.hdf_band_name:value(i-1)))
       end
    end
end

------------------------------------------------------------
--- Use radiance scaling sv fit but only if the ground type
--- is coxmunk
------------------------------------------------------------

ConfigCommon.radiance_scaling_sv_fit_coxmunk_only = ConfigCommon.radiance_scaling_sv_fit:new()
function ConfigCommon.radiance_scaling_sv_fit_coxmunk_only:create()
   local ground_type = self.config:ground_type_name()

   -- Only use this for coxmunk runs
   if(ground_type == "coxmunk") then
      self.config.using_radiance_scaling = true
      return ConfigCommon.radiance_scaling_sv_fit.create(self)
   else
      self.retrieve_bands = { false, false, false }
      return nil
   end
end

------------------------------------------------------------
--- Create a RadianceScalingLinearFit correction
------------------------------------------------------------

ConfigCommon.radiance_scaling_linear_fit = CreatorMultiSpec:new()

function ConfigCommon.radiance_scaling_linear_fit:create()
   local res = {}
   for i=1,self.config.number_pixel:rows() do
      local band_ref = self.config.common.band_reference
      local rad_meas = self.config.l1b:radiance_spectral_range(i-1)
      local fit_offset = self.fit_offset ~= nil and self.fit_offset[i]
      res[i] = RadianceScalingLinearFit(rad_meas,
                                        band_ref(i-1),
                                        self.config.common.desc_band_name:value(i-1),
                                        fit_offset)
   end
   self.radiance_scaling = res
   return res
end

function ConfigCommon.radiance_scaling_linear_fit:register_output(ro)
    if (self.radiance_scaling) then
       for i=1,self.config.number_pixel:rows() do
          ro:push_back(RadianceScalingOutput.create(self.radiance_scaling[i], 
                       self.config.common.hdf_band_name:value(i-1)))
       end
    end
end

function ConfigCommon.radiance_scaling_linear_fit:initial_guess(i)
   -- This class never uses the state vector
   return nil 
end

------------------------------------------------------------
--- Use radiance scaling linear fit but only if the ground type
--- is coxmunk
------------------------------------------------------------

ConfigCommon.radiance_scaling_linear_fit_coxmunk_only = ConfigCommon.radiance_scaling_linear_fit:new()
function ConfigCommon.radiance_scaling_linear_fit_coxmunk_only:create()
   local ground_type = self.config:ground_type_name()

   -- Only use this for coxmunk runs
   if(ground_type == "coxmunk") then
      self.config.using_radiance_scaling = true
      return ConfigCommon.radiance_scaling_linear_fit.create(self)
   else
      self.retrieve_bands = { false, false, false }
      return nil
   end
end

------------------------------------------------------------
--- Create a FluorescenceEffect effect
------------------------------------------------------------

ConfigCommon.fluorescence_effect = CreatorApriori:new()

function ConfigCommon.fluorescence_effect:create()
   local flag = self:retrieval_flag()

   local rad_unit = Unit(self.config.l1b:radiance_with_unit(0).units)
   self.fluorescence_effect = FluorescenceEffect(self:apriori(), flag,
                     self.config.atmosphere,
                     self.config.stokes_coefficient,
                     self.config.l1b:zen_with_unit(0),
                     0, self:reference_point(), rad_unit)
   local res = {}
   res[1] =  self.fluorescence_effect
   return res
end

function ConfigCommon.fluorescence_effect:initial_guess(i)
   if (i ~= 1 or not self.fluorescence_effect) then
       return nil
   end

   local rflag = self:retrieval_flag()
   local ig = InitialGuessValue()
   ig:apriori_subset(rflag, self:apriori_v())

   ig:apriori_covariance_subset(rflag, self:covariance())
   ig.initial_guess = self:iguess_v()
   return ig
end

function ConfigCommon.fluorescence_effect:register_output(ro)
    if(self.fluorescence_effect) then
        ro:push_back(FluorescenceEffectOutput.create(self.fluorescence_effect))
    end
end

------------------------------------------------------------
--- Apriori for fluorescence. We read the apriori from the
--- given HDF field. If we have the IMAP file as input, 
--- we change the first entry to use an initial guess based 
--- on the SIF results.
------------------------------------------------------------

function ConfigCommon.fluorescence_apriori(field)
   return function(self)
      local apr = self.config:h():apriori(field)
      if(self.config:h_imap()) then
	 local sid = self.config:l1b_sid_list()
	 local sounding_index = sid:sounding_number()
	 local idx = sid:frame_number()
	 local sif_757 = self.config:h_imap():read_double_2d("DOASFluorescence/fluorescence_radiance_757nm_corr_idp", idx, sounding_index)
	 local sif_771 = self.config:h_imap():read_double_2d("DOASFluorescence/fluorescence_radiance_771nm_corr_idp", idx, sounding_index)
	 if(sif_757 == -999999.0) then
	    -- if SIF757 failed, use only SIF771
	    apr:set(0, 1.8 * sif_771)
	 elseif(sif_771 == -999999.0) then
	    -- if SIF771 failed, use only SI757
	    apr:set(0, sif_757)
	 else
	    -- Average values
	    -- This value of 1.8 was supplied by Christian, see Ticket #2160
	    apr:set(0, (sif_757 + 1.8 * sif_771)/2)
	 end
      end
      return apr
   end
end

------------------------------------------------------------
--- Covariance for fluorescence. We read the covariance out
--- of the given HDF field. We then change the first entry
--- to:
---  1. If we have the IMAP file, we calculate the sigma of
---     the SIF and scale it
---  2. If we don't have IMAP, e take the HDF file and 
---     scale it by the radiance.
------------------------------------------------------------

function ConfigCommon.fluorescence_covariance(field)
   return function(self)
      local cov = self.config:h():covariance(field)
      if(self.config:h_imap()) then
	 local sid = self.config:l1b_sid_list()
	 local sounding_index = sid:sounding_number()
	 local idx = sid:frame_number()
	 local sif_757 = self.config:h_imap():read_double_2d("DOASFluorescence/fluorescence_radiance_757nm_corr_idp", idx, sounding_index)
	 local t = self.config:h_imap():read_double_2d("DOASFluorescence/fluorescence_radiance_757nm_uncert_idp", idx, sounding_index)
	 if(sif_757 == -999999.0) then
            -- if SIF757 failed, fall back to SIF771.
	    t = self.sif_uncert_ratio * self.config:h_imap():read_double_2d("DOASFluorescence/fluorescence_radiance_771nm_uncert_idp", idx, sounding_index)
	 end
	 if(t < 0) then
	    error("both IDP retrievals failed, cannot set fluorescence prior")
	 end
	 local sigma = self.sif_sigma_scale * t
	 cov:set(0, 0, sigma * sigma)
      else
	 -- Make Fs_755 covariance as as x% of the continuum level
	 local rad = self.config.l1b:radiance(0)
	 local nrad = rad:rows()
	 local scaled_cov = 
	    math.sqrt(cov(0,0)) * rad(Range(nrad-11,nrad-1)):mean()
	 cov:set(0, 0, scaled_cov * scaled_cov)
      end
      return cov
   end
end
------------------------------------------------------------
--- Fluorsecence only for lambertian runs
------------------------------------------------------------

ConfigCommon.fluorescence_effect_lambertian_only = ConfigCommon.fluorescence_effect:new()
function ConfigCommon.fluorescence_effect_lambertian_only:create()
   local ground_type = self.config:ground_type_name()

   -- Only use this for coxmunk runs
   if(ground_type == "lambertian") then
      return ConfigCommon.fluorescence_effect.create(self)
   else
      return nil
   end
end

------------------------------------------------------------
--- Use a precomputed noise file.
------------------------------------------------------------

ConfigCommon.noise_ascii = Creator:new()

function ConfigCommon.noise_ascii:create()
   local noise_file
   if(self.config.noise_file) then
      noise_file = HeritageFile(self.config.noise_file)
   else
      noise_file = HeritageFile(self.config.spectrum_file)
   end
   return PrecomputedNoiseModel(noise_file)
end

------------------------------------------------------------
--- Variation of noise_ascii that returns an array. This is
--- useful in testing, when we want to mix a HDF file wiht
--- an ASCII noise file.
------------------------------------------------------------

ConfigCommon.noise_ascii_array = Creator:new()

function ConfigCommon.noise_ascii_array:create()
   local noise_file
   if(self.config.noise_file) then
      noise_file = HeritageFile(self.config.noise_file)
   else
      noise_file = HeritageFile(self.config.spectrum_file)
   end
   local slist = self.config:l1b_sid_list()
   local res = {}
   for i=0,slist:size() - 1 do
      res[i] = PrecomputedNoiseModel(noise_file)
   end
   return res
end

------------------------------------------------------------
--- Creator for input objects, almost every mode has a
--- L1B file
------------------------------------------------------------

CreatorInput = CompositeCreator:new()

function CreatorInput:create_parent_object(sub_object)
    -- Noop
    return sub_object
end

function CreatorInput:sub_object_key()
   return {"l1b"}
end

ConfigCommon.l1b_input = CreatorInput:new()

ConfigCommon.l1b_met_input = CreatorInput:new()

function ConfigCommon.l1b_met_input:sub_object_key()
   return {"l1b", "met"}
end

------------------------------------------------------------
--- Common stuff in creating a L1b object
------------------------------------------------------------

CreatorL1b = CompositeCreator:new()

function CreatorL1b:sub_object_key()
   return {"noise"}
end

function CreatorL1b:register_output(ro)
   ro:push_back(Level1bOutput(self.config.l1b))
end

------------------------------------------------------------
--- Create a level 1b file were we read it from a ASCII file.
--- Real data comes from HDF, but we have old test data that
--- uses this format.
------------------------------------------------------------

ConfigCommon.level1b_ascii = CreatorL1b:new()

function ConfigCommon.level1b_ascii:create_parent_object()
   return Level1bHeritage(self.config.soundinginfo_file, self.config.spectrum_file,
                          self.config.noise)
end

------------------------------------------------------------
--- Creates a L1B object that scales another L1B objects'
--- radiance per spectral band
------------------------------------------------------------

ConfigCommon.level1b_scale_radiance = CreatorL1b:new()

function ConfigCommon.level1b_scale_radiance:create()
    -- Unset the register output of the unscaled L1B creator so it can
    -- be set by this one
    if (self.unscaled_l1b == nil) then
        error("unscaled_l1b block undefined for level1b_scale_radiance")
    end
    self.unscaled_l1b.creator.register_output = nil
    return CompositeCreator.create(self)
end

function ConfigCommon.level1b_scale_radiance:sub_object_key()
   return {"unscaled_l1b"}
end

function ConfigCommon.level1b_scale_radiance:create_parent_object()
    local nspec = self.config.spec_win:number_spectrometer()
    local scaling = Blitz_double_array_1d(nspec)
    for i, v in ipairs(self.scaling) do
        scaling:set(i-1, v)
    end
    return Level1bScaleRadiance(self.config.unscaled_l1b, scaling)
end

------------------------------------------------------------
--- Create a SolarAbsorptionAndContinuum solar model
------------------------------------------------------------

ConfigCommon.solar_absorption_and_continuum = CompositeCreator:new()

function ConfigCommon.solar_absorption_and_continuum:sub_object_key()
   return {"doppler_shift", "solar_absorption", "solar_continuum"}
end

function ConfigCommon.solar_absorption_and_continuum:create_parent_object(sub_object)
   local res = {}
   for i=1,self.config.number_pixel:rows() do
      res[i] = SolarAbsorptionAndContinuum(self.config.doppler_shift[i],
           self.config.solar_absorption[i], self.config.solar_continuum[i])
   end
   return res
end

------------------------------------------------------------
--- Create solar doppler shift from Level 1 data.
------------------------------------------------------------

ConfigCommon.solar_doppler_from_l1b = Creator:new()

function ConfigCommon.solar_doppler_from_l1b:create()
   local res = {}
   local i
   for i=1,self.config.number_pixel:rows() do
      if(self.config.l1b.has_solar_relative_velocity ~= nil and
	 self.config.l1b:has_solar_relative_velocity()) then
	 res[i] = SolarDopplerShiftL1b(self.config.l1b:solar_distance(),
			       self.config.l1b:solar_velocity(),
			       self.do_doppler_shift)
      else
	 res[i] = SolarDopplerShiftPolynomial.create_from_l1b(self.config.l1b, 
           i - 1, self.do_doppler_shift)
      end
   end
   return res
end

------------------------------------------------------------
--- Create solar doppler shift from FTS runlog.
------------------------------------------------------------

ConfigCommon.solar_doppler_from_runlog = Creator:new()

function ConfigCommon.solar_doppler_from_runlog:create()
   local res = {}
   local i
   for i=1,self.config.number_pixel:rows() do
      res[i] = SolarDopplerShiftPolynomial.create_from_runlog(self.config.l1b, 
           i - 1, self.do_doppler_shift)
   end
   return res
end

------------------------------------------------------------
--- Create SolarAbsorptionTable
------------------------------------------------------------

ConfigCommon.solar_absorption_table = Creator:new()

function ConfigCommon.solar_absorption_table:create()
   local res = {}
   local i
   for i=1,self.config.number_pixel:rows() do
      res[i] = SolarAbsorptionTable(self.config:h_solar(),
                                    "/Solar/Absorption/Absorption_" .. i)
   end
   return res
end

--- This is replaced with SolarAbsorptionTable. We'll leave the code
--- here for a little while, in case we come across some odd requirement
--- to be able to use this old code again (e.g., compare with an old run)

--- ConfigCommon.solar_absorption_oco_file = Creator:new()

--- function ConfigCommon.solar_absorption_oco_file:create()
---   local res = {}
---   local i
---   for i=1,self.config.number_pixel:rows() do
---      res[i] = SolarAbsorptionOcoFile(self.config:h(), self.hdf_group)
---   end
---   return res
--- end

------------------------------------------------------------
--- Create SolarContinuumTable
------------------------------------------------------------

ConfigCommon.solar_continuum_table = Creator:new()

function ConfigCommon.solar_continuum_table:create()
   local res = {}
   local i
   for i=1,self.config.number_pixel:rows() do
      res[i] = SolarContinuumTable(self.config:h_solar(),
                                   "/Solar/Continuum/Continuum_" .. i,
                                   self.convert_from_photon)
   end
   return res
end

------------------------------------------------------------
--- Create SolarContinuumPolynomial.
------------------------------------------------------------

--- This is replaced with SolarContinuumTable. We'll leave the code
--- here for a little while, in case we come across some odd requirement
--- to be able to use this old code again (e.g., compare with an old run)

-- ConfigCommon.solar_continuum_polynomial = Creator:new()

-- function ConfigCommon.solar_continuum_polynomial:create()
--   local res = {}
--   local i
--   for i=1,self.config.number_pixel:rows() do
--      res[i] = SolarContinuumPolynomial(self:apriori(), 
--                                        self.convert_from_photon)
--   end
--   return res
-- end

------------------------------------------------------------
--- Creates an instrument doppler spectrum effect that uses
--- the relative velocity from the L1B file to shift the
--- high resolution spectral grid accordingly.
------------------------------------------------------------

ConfigCommon.instrument_doppler = CreatorApriori:new()

function ConfigCommon.instrument_doppler:create()
    local res = {}
    local rel_vel = self:apriori_v()
    local flag = self:retrieval_flag()
    for i=1,self.config.number_pixel:rows() do
       res[i] = InstrumentDoppler(rel_vel(i-1), "m / s", flag(i-1)) 
    end
    return res
end

function ConfigCommon.instrument_doppler:apriori_v()
    local nspec = self.config.number_pixel:rows()
    local ap = Blitz_double_array_1d(nspec)
    for i=1,nspec do
        ap:set(i-1, self.config.l1b:relative_velocity(i-1))
    end
    return ap
end

function ConfigCommon.instrument_doppler:covariance_v()
    if self.covariance ~= nil then
        return self:covariance()
    else
        local nspec = self.config.number_pixel:rows()
        local cov = Blitz_double_array_2d(nspec,nspec)
        for i=1,nspec do
            cov:set(i-1, i-1, 1e-20)
        end
        return cov
    end
end

------------------------------------------------------------
--- Create overall atmosphere. 
------------------------------------------------------------

ConfigCommon.atmosphere_oco = CompositeCreator:new()

function ConfigCommon.atmosphere_oco:sub_object_key()
-- Order here is important, because some objects require other objects
-- be created first (e.g., most everything requires "pressure"
   return {"constants", "pressure", "temperature", "ground", 
           "altitude", "absorber", "relative_humidity", "aerosol"}
end

function ConfigCommon.atmosphere_oco:sub_initial_guess_key()
-- Order stuff goes in the state vector. This is arbitrary, but this
-- is the "expected" order because it was what was done historically. 
-- This isn't the same order as we create objects, again just because
-- of history.
   return {"absorber", "pressure", "temperature", "aerosol", "ground"}
end

function ConfigCommon.atmosphere_oco:create_parent_object(sub_object)
   local c = self.config

   c.number_pressure_level = c.pressure:max_number_level()
   if(c.aerosol and c.ground) then
      return AtmosphereOco(c.absorber, c.pressure, c.temperature, c.aerosol, 
                           c.relative_humidity, c.ground, c.altitude, 
			   c.constants)
   elseif(c.aerosol) then
      return AtmosphereOco(c.absorber, c.pressure, c.temperature, c.aerosol, 
                           c.relative_humidity, c.altitude, c.constants)
   elseif(c.ground) then
      return AtmosphereOco(c.absorber, c.pressure, c.temperature,  
                           c.relative_humidity, c.ground, c.altitude, 
			   c.constants)
   end
   return AtmosphereOco(c.absorber, c.pressure, c.temperature, 
			c.relative_humidity, c.altitude,
                        c.constants)
end

function ConfigCommon.atmosphere_oco:add_to_statevector(sv)
   CompositeCreator.add_to_statevector(self, sv)
   sv:add_observer(self.config.atmosphere)
end

------------------------------------------------------------
--- Create pressure with the fixed levels, retrieving surface
--- pressure and initial guess.
--- 
--- This also has the side effect of creating config.pinp which
--- can then be used by other classes.
------------------------------------------------------------

ConfigCommon.pressure_fixed_level = CreatorApriori:new {}

function ConfigCommon.pressure_fixed_level:create()
   self.config.pinp = PressureLevelInput(self:pressure_levels())
   return PressureFixedLevel(self:retrieval_flag()(0), self.config.pinp, 
                             self:apriori()(0))
end

function ConfigCommon.pressure_fixed_level:register_output(ro)
   ro:push_back(PressureFixedLevelOutput.create(self.config.pressure,
                                                self.config.state_vector))
end

------------------------------------------------------------
--- Create pressure with the sigma levels, retrieving surface
--- pressure and initial guess.
------------------------------------------------------------

ConfigCommon.pressure_sigma = CreatorApriori:new {}

function ConfigCommon.pressure_sigma:create()
   return PressureSigma(self:a(), self:b(), self:apriori()(0),
                        self:retrieval_flag()(0))
end

function ConfigCommon.pressure_sigma:register_output(ro)
   ro:push_back(PressureOutput(self.config.pressure,
                               self.config.state_vector))
end

------------------------------------------------------------
--- Create pressure with the sigma levels, determining
--- the sigma levels a and b from an input pressure grid.
--- Retrieves surface pressure
------------------------------------------------------------

ConfigCommon.pressure_sigma_profile = CreatorApriori:new {}

function ConfigCommon.pressure_sigma_profile:create()
   return PressureSigma(self:pressure_levels(), self:apriori()(0),
                        self:retrieval_flag()(0))
end

function ConfigCommon.pressure_sigma_profile:register_output(ro)
   ro:push_back(PressureOutput(self.config.pressure,
                            self.config.state_vector))
end

------------------------------------------------------------
--- Temperature at fixed levels, where we just fit for the
--- offset.
------------------------------------------------------------

ConfigCommon.temperature_fixed_level = Creator:new {}

function ConfigCommon.temperature_fixed_level:retrieval_flag()
   flag = Blitz_bool_array_1d(1)

   if self.retrieved ~= nil then
      flag:set(Range.all(), self.retrieved)
   else
      flag:set(Range.all(), true)
   end

   return flag
end

function ConfigCommon.temperature_fixed_level:create()
   self.levels.config = self.config
   self.offset.config = self.config
   local tlev = self.levels:apriori()
   local toff = self.offset:apriori()(0)
   local flag_temp = Blitz_bool_array_1d(tlev:rows())
   flag_temp:set(Range.all(), false)
   local flag_offset = self:retrieval_flag()(0)
   return TemperatureFixedLevel(flag_temp, flag_offset, tlev, toff, 
                                self.config.pressure, self.config.pinp)
end

function ConfigCommon.temperature_fixed_level:initial_guess()
   local ig = InitialGuessValue()
   self.offset.config = self.config
   ig:apriori_subset(self:retrieval_flag(), self.offset:apriori())
   ig:apriori_covariance_subset(self:retrieval_flag(), self.offset:covariance())

   return ig
end

function ConfigCommon.temperature_fixed_level:register_output(ro)
   ro:push_back(TemperatureFixedLevelOutput.create(self.config.temperature))
end

------------------------------------------------------------
--- Temperature using ECMWF, where we just fit for the
--- offset.
------------------------------------------------------------

ConfigCommon.temperature_met = CreatorApriori:new {}

function ConfigCommon.temperature_met:create()
   return TemperatureMet(self.config.met, self.config.pressure,
                         self:apriori()(0), self:retrieval_flag()(0))
end

function ConfigCommon.temperature_met:register_output(ro)
   ro:push_back(TemperatureMetOutput.create(self.config.temperature))
end

------------------------------------------------------------
--- Temperature using specificed level values
--- where we just fit for the offset.
------------------------------------------------------------

ConfigCommon.temperature_level_offset = CreatorApriori:new {}

function ConfigCommon.temperature_level_offset:create()
   return TemperatureLevelOffset(self.config.pressure, self:temperature_levels(),
                                 self:apriori()(0), self:retrieval_flag()(0))
end

function ConfigCommon.temperature_level_offset:register_output(ro)
   ro:push_back(TemperatureLevelOffsetOutput.create(self.config.temperature))
end

------------------------------------------------------------
--- Temperature using specificed level values
--- where we fit all levels
------------------------------------------------------------

ConfigCommon.temperature_level = CreatorApriori:new {}

function ConfigCommon.temperature_level:create()
   return TemperatureLevel(self:apriori(), self.config.pressure, self:retrieval_flag())
end

function ConfigCommon.temperature_level:register_output(ro)
  --ro:push_back(TemperatureLevelOutput.create(self.config.temperature))
end

------------------------------------------------------------
--- Lambertian ground state vector component and initial guess
------------------------------------------------------------

ConfigCommon.lambertian_retrieval = CreatorMultiSpec:new {}

function ConfigCommon.lambertian_retrieval:create()
   local num_coeff = self:apriori_v(0):rows()
   local num_spec = self.config.number_pixel:rows()

   local ap = Blitz_double_array_2d(num_spec, num_coeff)
   local flag = Blitz_bool_array_2d(num_spec, num_coeff)

   for i = 1, num_spec do
       ap:set(i-1, Range.all(), self:apriori_v(i - 1))
       flag:set(i-1, Range.all(), self:retrieval_flag(i))
   end

   local lambertian = GroundLambertian(ap, flag, 
                                       self.config.common.band_reference,
                                       self.config.common.desc_band_name)
   return lambertian
end

function ConfigCommon.lambertian_retrieval:initial_guess()
   local ig = CompositeInitialGuess()
   for i=1,self.config.number_pixel:rows() do
       local flag = self:retrieval_flag(i)

       local ap = self:apriori_v(i - 1) 
       local apcov_in = self:covariance_v(i - 1)

       -- Copy input covariance to modified one which repeats the last covariance value if
       -- the number of degrees is defined
       local apcov_mod = Blitz_double_array_2d(ap:rows(), ap:rows())
       apcov_mod:set(Range.all(), Range.all(), 0.0)
       src_idx = 0
       for dst_idx = 0, ap:rows()-1 do
          apcov_mod:set(dst_idx, dst_idx, apcov_in(src_idx, src_idx))
          if dst_idx < apcov_in:rows()-1 then
             src_idx = src_idx + 1
          end
       end

       local band_ig = InitialGuessValue()
       band_ig:apriori_subset(flag, ap)
       band_ig:apriori_covariance_subset(flag, apcov_mod)
       ig:add_builder(band_ig)
   end
   return ig 
end

------------------------------------------------------------
--- Lambertian ground state vector component and initial guess
------------------------------------------------------------

ConfigCommon.brdf_scale_retrieval = CreatorMultiSpec:new {}

function ConfigCommon.brdf_scale_retrieval:create()
   local num_coeff = self:apriori_v(0):rows()
   local num_spec = self.config.number_pixel:rows()

   local ap = Blitz_double_array_2d(num_spec, num_coeff)
   local flag = Blitz_bool_array_2d(num_spec, num_coeff)

   for i = 1, num_spec do
       ap:set(i-1, Range.all(), self:apriori_v(i - 1))
       flag:set(i-1, Range.all(), self:retrieval_flag(i))
   end

   local lambertian = GroundBrdfWeight(ap, flag, 
                                       self.config.common.band_reference,
                                       self.config.common.desc_band_name,
                                       self.scaled_brdf_name)
   return lambertian
end

function ConfigCommon.brdf_scale_retrieval:initial_guess()
   local ig = CompositeInitialGuess()
   for i=1,self.config.number_pixel:rows() do
       local flag = self:retrieval_flag(i)

       local ap = self:apriori_v(i - 1) 
       local apcov_in = self:covariance_v(i - 1)

       -- Copy input covariance to modified one which repeats the last covariance value if
       -- the number of degrees is defined
       local apcov_mod = Blitz_double_array_2d(ap:rows(), ap:rows())
       apcov_mod:set(Range.all(), Range.all(), 0.0)
       src_idx = 0
       for dst_idx = 0, ap:rows()-1 do
          apcov_mod:set(dst_idx, dst_idx, apcov_in(src_idx, src_idx))
          if dst_idx < apcov_in:rows()-1 then
             src_idx = src_idx + 1
          end
       end

       local band_ig = InitialGuessValue()
       band_ig:apriori_subset(flag, ap)
       band_ig:apriori_covariance_subset(flag, apcov_mod)
       ig:add_builder(band_ig)
   end
   return ig 
end

------------------------------------------------------------
--- Creates a GroundOco only using a Lambertian retrieval
------------------------------------------------------------

ConfigCommon.ground_lambertian = CompositeCreator:new {}

function ConfigCommon.ground_lambertian:sub_object_key()
   return { "lambertian" }
end

function ConfigCommon.ground_lambertian:create_parent_object(sub_object)
   return self.config.lambertian
end

function ConfigCommon.ground_lambertian:register_output(ro)
   ro:push_back(GroundLambertianOutput.create(self.config.lambertian, self.config.common.hdf_band_name))
end

------------------------------------------------------------
--- Coxmunk ground state vector component and initial guess
------------------------------------------------------------

ConfigCommon.coxmunk_retrieval = CreatorApriori:new {}

function ConfigCommon.coxmunk_retrieval:create()
   return GroundCoxmunk(self:apriori_v()(0), self:retrieval_flag()(0),
                        self:refractive_index())
end

------------------------------------------------------------
--- Creates a GroundOco only using a Coxmunk retrieval
------------------------------------------------------------

ConfigCommon.ground_coxmunk = CompositeCreator:new {}

function ConfigCommon.ground_coxmunk:sub_object_key()
   return { "coxmunk" }
end

function ConfigCommon.ground_coxmunk:create_parent_object(sub_object)
   return self.config.coxmunk
end

function ConfigCommon.ground_coxmunk:register_output(ro)
   ro:push_back(GroundCoxmunkOutput.create(self.config.coxmunk))
end

------------------------------------------------------------
--- Creates a GroundOco using both Lambertian and Coxmunk 
--- retrieval
------------------------------------------------------------

ConfigCommon.ground_coxmunk_plus_lamb = CompositeCreator:new {}

function ConfigCommon.ground_coxmunk_plus_lamb:sub_object_key()
   -- Ordering is important here or else the statevector is setup incorrectly
   return { "coxmunk", "coxmunk_lambertian" }
end

function ConfigCommon.ground_coxmunk_plus_lamb:create_parent_object(sub_object)
   return GroundCoxmunkPlusLambertian.create(self.config.coxmunk, self.config.coxmunk_lambertian)
end

function ConfigCommon.ground_coxmunk_plus_lamb:register_output(ro)
   ro:push_back(GroundCoxmunkPlusLambertianOutput.create(self.config.coxmunk, self.config.coxmunk_lambertian, self.config.common.hdf_band_name))
end

------------------------------------------------------------
--- Creates a GroundOco using both Lambertian and Coxmunk 
--- retrieval
------------------------------------------------------------

ConfigCommon.ground_coxmunk_scaled = CompositeCreator:new {}

function ConfigCommon.ground_coxmunk_scaled:sub_object_key()
   -- Ordering is important here or else the statevector is setup incorrectly
   return { "coxmunk", "coxmunk_scaled" }
end

function ConfigCommon.ground_coxmunk_scaled:create_parent_object(sub_object)
   return GroundCoxmunkScaled.create(self.config.coxmunk, self.config.coxmunk_scaled)
end

function ConfigCommon.ground_coxmunk_scaled:register_output(ro)
   ro:push_back(GroundCoxmunkScaledOutput.create(self.config.l1b, self.config.coxmunk, self.config.coxmunk_scaled, self.config.common.hdf_band_name))
end

------------------------------------------------------------
--- Modify the GroundBrdf a_priori, scaling by the 
--- kernel value
------------------------------------------------------------

function ConfigCommon.brdf_weight(self, brdf_class, ap, i)
   local sza_d = self.config.l1b:sza()(i) 
   local vza_d = self.config.l1b:zen()(i) 
   local azm_d = self.config.l1b:azm()(i) 
   local alb_cont = self.config:albedo_from_signal_level(1)(self, i)(0)

   -- Extract all but the slope portion of the apriori to feed into the
   -- albedo calculation function
   local params = Blitz_double_array_1d(5)
   params:set(Range.all(), ap(Range(0, 4)))

   local alb_calc = brdf_class.kernel_value(params, sza_d, vza_d, azm_d)
   local weight = alb_cont / alb_calc

   return weight
end

function ConfigCommon.brdf_veg_apriori(field)
    return function(self, i)
        local ap = self.config:h():apriori(field, i) 
        local weight = ConfigCommon.brdf_weight(self, GroundBrdfVeg, ap, i)
        ap:set(5, ap(5) * weight)
        return ap
    end
end

function ConfigCommon.brdf_soil_apriori(field)
    return function(self, i)
        local ap = self.config:h():apriori(field, i) 
        local weight = ConfigCommon.brdf_weight(self, GroundBrdfSoil, ap, i)
        ap:set(5, ap(5) * weight)
        return ap
    end
end

------------------------------------------------------------
--- Common behavior of both BRDF retrieval methods       ---
------------------------------------------------------------

ConfigCommon.brdf_retrieval = CreatorMultiSpec:new {}

function ConfigCommon.brdf_retrieval:retrieval_flag(i)
   local flag = Blitz_bool_array_1d(self:apriori_v(i - 1):rows())

   n_coefs = self:apriori_v(0):rows()

   if self.retrieve_bands ~= nil and self.retrieve_bands[i] then
       flag:set(Range.all(), false)
       for i = 5, n_coefs - 1 do
           flag:set(i, true)
       end
   else
        flag:set(Range.all(), false)
   end

   return flag
end

------------------------------------------------------------
--- Breon veg ground state vector component and initial guess
------------------------------------------------------------

ConfigCommon.brdf_veg_retrieval = ConfigCommon.brdf_retrieval:new {}

function ConfigCommon.brdf_veg_retrieval:create()
   local num_spec = self.config.number_pixel:rows()

   n_coefs = self:apriori_v(0):rows()

   local ap = Blitz_double_array_2d(num_spec, n_coefs)
   local flag = Blitz_bool_array_2d(num_spec, n_coefs)

   for i = 1, num_spec do
       ap:set(i-1, Range.all(), self:apriori_v(i - 1))
       flag:set(i-1, Range.all(), self:retrieval_flag(i))
   end

   return GroundBrdfVeg(ap, flag, self.config.common.band_reference, self.config.common.desc_band_name)
end

------------------------------------------------------------
--- Sets up the ground using only using a Breon veg retrieval
------------------------------------------------------------

ConfigCommon.ground_brdf_veg = CompositeCreator:new {}

function ConfigCommon.ground_brdf_veg:sub_object_key()
   return { "brdf_veg" }
end

function ConfigCommon.ground_brdf_veg:create_parent_object(sub_object)
   return self.config.brdf_veg
end

function ConfigCommon.ground_brdf_veg:register_output(ro)
   ro:push_back(GroundBrdfOutput.create(self.config.brdf_veg, self.config.l1b, self.config.common.hdf_band_name))
end

------------------------------------------------------------
--- Breon soil ground state vector component and initial guess
------------------------------------------------------------

ConfigCommon.brdf_soil_retrieval = ConfigCommon.brdf_retrieval:new {}

function ConfigCommon.brdf_soil_retrieval:create()
   local num_spec = self.config.number_pixel:rows()

   n_coefs = self:apriori_v(0):rows()

   local ap = Blitz_double_array_2d(num_spec, n_coefs)
   local flag = Blitz_bool_array_2d(num_spec, n_coefs)

   for i = 1, num_spec do
       ap:set(i-1, Range.all(), self:apriori_v(i - 1))
       flag:set(i-1, Range.all(), self:retrieval_flag(i))
   end

   return GroundBrdfSoil(ap, flag, self.config.common.band_reference, self.config.common.desc_band_name)
end

------------------------------------------------------------
--- Sets up the ground using only using a Breon soil retrieval
------------------------------------------------------------

ConfigCommon.ground_brdf_soil = CompositeCreator:new {}

function ConfigCommon.ground_brdf_soil:sub_object_key()
   return { "brdf_soil" }
end

function ConfigCommon.ground_brdf_soil:create_parent_object(sub_object)
   return self.config.brdf_soil
end

function ConfigCommon.ground_brdf_soil:register_output(ro)
   ro:push_back(GroundBrdfOutput.create(self.config.brdf_soil, self.config.l1b, self.config.common.hdf_band_name))
end

------------------------------------------------------------
--- Create lambertian ground initial guess from radiance
--- and the other parts from the static HDF file
---
--- Requires arguments specific to the instrument the
--- radiances come from.
--- 
--- solar_strength - Strength of sun in radiance units at the band 
--- band_coninuums - Pixel ranges where the continuum radiance
---   can be found
--- use_range_max - true or false per spectrometer on
---   whether to use the maximum value of the continuum
---   range instead of including all points in the average
---
--- Example from GOSAT:
---      SOLAR_STRENGTH = {7.2e-6, 6.5e-6, 4.5e-6}
---      CONTINUUM_POINTS = 
---         { { Range(490, 529), Range(1586, 1602) },
---           { Range(2390, 2394), Range(2400, 2404) },
---           { Range(333, 336), Range(518, 521), Range(729, 732) },
---         }
---      USE_RANGE_MAX = { false, false, true }
------------------------------------------------------------

function ConfigCommon.calculate_albedo_from_radiance_constant(radiance, sza_r, solar_strength, continuum_points, use_range_max)
   local mean_count = 0
   local mean_sum = 0
   for idx, cont_range in ipairs(continuum_points) do
      -- Make sure we always get back an array object
      if type(cont_range) == "number" then
         cont_data = radiance(Range(cont_range,cont_range))
      else
         cont_data = radiance(cont_range)
      end

      if use_range_max then
         mean_sum = mean_sum + cont_data:max()
         mean_count = mean_count + 1
      else
         mean_sum = mean_sum + cont_data:sum()
         mean_count = mean_count + cont_data:rows()
      end
   end
   local albedo_val = math.pi * (mean_sum / mean_count) / (math.cos(sza_r) * solar_strength)

   return albedo_val
end

function ConfigCommon.calculate_albedo_from_radiance_polynomial(radiance, sza_r, solar_strength, continuum_points, use_range_max, polynomial_degree)
   local offset
   if (type(continuum_points[1]) == "table") then
      -- Case where seperate ranges for either end of the band are specified
      beg_alb = ConfigCommon.calculate_albedo_from_radiance_constant(radiance, sza_r, solar_strength, continuum_points[1], use_range_max)
      end_alb = ConfigCommon.calculate_albedo_from_radiance_constant(radiance, sza_r, solar_strength, continuum_points[2], use_range_max)
      offset = (beg_alb + end_alb) / 2.0
   else
      -- Just one set of points
      offset = ConfigCommon.calculate_albedo_from_radiance_constant(radiance, sza_r, solar_strength, continuum_points, use_range_max)
   end

   if not polynomial_degree then
      polynomial_degree = 1
   end

   local albedo_val = Blitz_double_array_1d(polynomial_degree + 1)
   albedo_val:set(0, offset)

   -- Set up albedo apriori values with proper number of polynomial degrees
   for d_idx = 1, polynomial_degree do
       albedo_val:set(d_idx, 0.0)
   end

   return albedo_val
end

function ConfigCommon:albedo_from_radiance(band_idx, solar_strength, continuum_points, use_range_max, polynomial_degree)
   if not solar_strength or not continuum_points or not use_range_max then
      error("Not all necessary arguments supplied")
   end

   l1b = self.config.l1b
   sza_d = l1b:sza()

   if not polynomial_degree then
      polynomial_degree = 1
   end

   sza_r = sza_d(band_idx) * math.pi / 180.0
   radiance = l1b:radiance(band_idx)

   albedo_val = ConfigCommon.calculate_albedo_from_radiance_polynomial(radiance, sza_r, 
                                                                       solar_strength[band_idx+1], 
                                                                       continuum_points[band_idx+1], 
                                                                       use_range_max[band_idx+1],
                                                                       polynomial_degree)
   return albedo_val
end

------------------------------------------------------------
-- Continuum signal level of measurement spectral 
-- Will only take into account samples used in retrieval
-- when called after spec_win, dispesion and l1b are set up
-- Indexing is zero based
------------------------------------------------------------

function ConfigCommon:meas_cont_signal(spec_idx)
    local signal
    if self.dispersion[spec_idx+1] ~= nil then
        local pixel_grid = self.dispersion[spec_idx+1]:pixel_grid()
        local grid_indexes = self.spec_win:grid_indexes(pixel_grid, spec_idx)
        signal = self.l1b:signal(spec_idx, grid_indexes) 
    else
        self:diagnostic_message("Not removing bad samples or out of range samples for measurement continuum signal calculation")
        signal = self.l1b:signal(spec_idx) 
    end
    return signal
end

------------------------------------------------------------
--- Create lambertian ground initial guess from radiance
--- using the signal level reported by the L1b object
--- Meant to replace albedo_from_radiance() above.
------------------------------------------------------------

function ConfigCommon:albedo_from_signal_level(polynomial_degree)
    if not polynomial_degree then
        polynomial_degree = 1
    end

    return function(self, spec_idx)
        local signal = self.config:meas_cont_signal(spec_idx).value
        local solar_strength = self.config.fm.atmosphere.ground.solar_strength[spec_idx+1]
        local sza_r = self.config.l1b:sza()(spec_idx) * math.pi / 180.0

        -- Account for solar distance Fsun = Fsun0 / (solar_distance_meters/AU)^2
        -- Create SolarDopplerShiftPolynomial so we can compute solar distance
        local solar_doppler_shift = SolarDopplerShiftPolynomial.create_from_l1b(self.config.l1b, spec_idx, true)
        local solar_dist = solar_doppler_shift:solar_distance().value
        solar_strength = solar_strength / solar_dist^2
     
        -- Account for stokes element for I
        local stokes_coef = self.config.l1b:stokes_coef()
        solar_strength = solar_strength * stokes_coef(spec_idx, 0)

        local offset = math.pi * signal / (math.cos(sza_r) * solar_strength)

        local albedo_val = Blitz_double_array_1d(polynomial_degree + 1)
        albedo_val:set(Range.all(), 0)
        albedo_val:set(0, offset)

        return albedo_val
    end
end

------------------------------------------------------------
--- Use a fixed albedo as a initial guess. Useful when 
--- matching simulated data where we know what the albedo
--- should be.
--- Example:
---   apriori = ConfigCommon.fixed_albedo({0.3, 0.3, 0.3}, {0.0, 0.0, 0.0}),
------------------------------------------------------------

function ConfigCommon.fixed_albedo(val1, val2)
   return function(self, i)
      local ap = Blitz_double_array_1d(2)
      ap:set(0,val1[i+1])
      ap:set(1,val2[i+1])
      return ap
   end
end

------------------------------------------------------------
--- For instances where a ground is not needed
------------------------------------------------------------

ConfigCommon.no_ground = Creator:new {}

function ConfigCommon.no_ground:create()
   return nil
end

------------------------------------------------------------
--- Creation for a aerosol is a little more complicated.
--- This base class handles that complication, derived classes
--- just need to supply the aerosol extinction.
------------------------------------------------------------

CreatorAerosol = CreatorApriori:new { }

function CreatorAerosol:apriori_v()
   return self:apriori(self.name)
end

function CreatorAerosol:covariance_v()
   return self:covariance(self.name)
end

function CreatorAerosol:create()
   local res = {}
   res.property = self:property(self.name)
   res.extinction = self:extinction()
   res.initial_guess = self:initial_guess()
   return res
end

------------------------------------------------------------
--- Create aerosol extinction model using logarithmic
--- profiles of exitinction.
------------------------------------------------------------

ConfigCommon.aerosol_log_profile = CreatorAerosol:new()

function ConfigCommon.aerosol_log_profile:extinction()
   return AerosolExtinctionLog(self.config.pressure, self:retrieval_flag(), 
                               self:apriori_v(), self.name)
end

------------------------------------------------------------
--- Create aerosol extinction model using linear
--- profiles of exitinction.
------------------------------------------------------------

ConfigCommon.aerosol_linear_profile = CreatorAerosol:new()

function ConfigCommon.aerosol_linear_profile:extinction()
   return AerosolExtinctionLinear(self.config.pressure, self:retrieval_flag(), 
                                  self:apriori_v(), self.name)
end

------------------------------------------------------------
--- Create aerosol extinction model using logarithmic
--- gaussian shape
------------------------------------------------------------

ConfigCommon.aerosol_log_shape_gaussian = CreatorAerosol:new()

function ConfigCommon.aerosol_log_shape_gaussian:extinction()
   return AerosolShapeGaussian(self.config.pressure, self:retrieval_flag(), 
                               self:apriori_v(), self.name, false)
end

------------------------------------------------------------
--- Create aerosol extinction model using linear
--- gaussian shape
------------------------------------------------------------

ConfigCommon.aerosol_linear_shape_gaussian = CreatorAerosol:new()

function ConfigCommon.aerosol_linear_shape_gaussian:extinction()
   return AerosolShapeGaussian(self.config.pressure, self:retrieval_flag(), 
                               self:apriori_v(), self.name, true)
end

------------------------------------------------------------
--- Create aerosol.
------------------------------------------------------------

ConfigCommon.aerosol_creator = CompositeCreator:new()

function ConfigCommon.aerosol_creator:sub_object_key()
   return self.aerosols
end

function ConfigCommon.aerosol_creator:create_parent_object(sub_object)
   self.vap = VectorAerosolProperty()
   self.vex = VectorAerosolExtinction()
   local i, t
   self.config.number_aerosol = 0
   for i, t in ipairs(sub_object) do
      self.vap:push_back(t.property)
      self.vex:push_back(t.extinction)
      self.config.number_aerosol = self.config.number_aerosol + 1
   end
   return AerosolOptical(self.vex, self.vap, self.config.pressure, 
			 self.config.relative_humidity)
end

function ConfigCommon.aerosol_creator:register_output(ro)
   ro:push_back(AerosolAodOutput(self.config.aerosol))

   for i=0, self.vex:size() - 1 do
      ro:push_back(AerosolParamOutput.create(self.vex:value(i)))
   end
end

------------------------------------------------------------
--- Create aerosol using the MerraAerosol class.
------------------------------------------------------------

ConfigCommon.merra_aerosol_creator = CompositeCreator:new()

function ConfigCommon:merra_file()
   if(os.getenv("merradir") and not self.merra_file_v) then
      local fname = string.format("%s/MERRA_Composite_Selection_M%02d_O2A.hdf5", 
				  os.getenv("merradir"), self.l1b:month(0))
      self.merra_file_v = HdfFile(fname)
      self.input_file_description = self.input_file_description .. 
	 "Merra input file:    " .. fname .. "\n"
      self.merra_file_name = fname
   else if(self.merra_dir and not self.merra_file_v) then
	 local fname = string.format("%s/MERRA_Composite_Selection_M%02d_O2A.hdf5", 
				     self.merra_dir, self.l1b:month(0))
	 self.merra_file_v = HdfFile(fname)
	 self.input_file_description = self.input_file_description .. 
	    "Merra input file:    " .. fname .. "\n"
	 self.merra_file_name = fname
      end
   end
   if(self.merra_file_v) then
      return self.merra_file_v
   end
   error({code=-1})
end

function ConfigCommon.merra_aerosol_creator:sub_object_key()
   return self.aerosols
end

function ConfigCommon.merra_aerosol_creator:create_parent_object(sub_object)
   -- Luabind can only handle up to 10 arguments per function. As an easy
   -- work around we put the various thresholds into an array.
   local mq = Blitz_double_array_1d(7)
   mq:set(0, self.max_aod)
   mq:set(1, self.exp_aod)
   mq:set(2, self.min_types)
   mq:set(3, self.max_types)
   if(self.linear_aod) then
      mq:set(4, 1)
   else
      mq:set(4, 0)
   end
   if(self.relative_humidity_aerosol) then
      mq:set(5, 1)
   else
      mq:set(5, 0)
   end
   mq:set(6, self.max_residual)

   self.merra_aerosol = MerraAerosol.create(
     self.config:merra_file(), self.config:h_merra_aerosol(),
     self.config.l1b:latitude(0),
     self.config.l1b:longitude(0),
     self.config.pressure,
     self.config.relative_humidity,
     -- In production, have a single covariance for all merra types,
     -- but allow for using functions that have a different covariance
     -- for different types
     self:covariance(0), 
     mq)
   self.sub_object_save = sub_object
   local i, t
   for i, t in ipairs(sub_object) do
      self.merra_aerosol:add_aerosol(t.extinction, t.property, t.initial_guess)
   end
   local res = self.merra_aerosol:aerosol()
   self.config.number_aerosol = res:number_particle()
   return res
end

function ConfigCommon.merra_aerosol_creator:initial_guess()
   local cig = CompositeInitialGuess()
   local flag = ConfigCommon.create_flag(3, Range.all())
   local cig_merra = CompositeInitialGuess()
   cig_merra:add_builder(self.merra_aerosol:initial_guess())

   --- Use AOD from merra, but leave the remaining shape parameters at values
   --- fixed at the supplied apriori values (e.g., read from HDF file, or 
   --- whatever).
   local ig_merra = cig_merra:initial_guess()
   for i=1,self.merra_aerosol:number_merra_particle() do
      -- Send particle number to ap and cov functions in case they
      -- can use them
      local oval = self:apriori(i-1)

      -- Optionally ignore the merra initial guess value and just use
      -- value from apriori function
      if (self.ignore_merra_aod == nil) then
        oval:set(0, ig_merra((i - 1) * 3))
      end

      local ig = InitialGuessValue()
      ig:apriori_subset(flag, oval)
      ig:apriori_covariance_subset(flag, self:covariance(i-1))
      cig:add_builder(ig)
   end

   --- Now add in the fixed particles.
   local i, t
   for i, t in ipairs(self.sub_object_save) do
      cig:add_builder(t.initial_guess)
   end

   return cig
end

function ConfigCommon.merra_aerosol_creator:register_output(ro)
   local all_aer_names = self.config:merra_file():read_string_vector("/COMPOSITE_NAME")
   for i, aer_name in ipairs(self.aerosols) do
       all_aer_names:push_back(aer_name)
   end
 
   ro:push_back(AerosolConsolidatedOutput(self.config.aerosol, all_aer_names))
end

------------------------------------------------------------
--- Create aerosol using the AerosolMetPrior class.
------------------------------------------------------------

ConfigCommon.aerosol_met_prior_creator = CompositeCreator:new()

function ConfigCommon.aerosol_met_prior_creator:sub_object_key()
   return self.aerosols
end

function ConfigCommon.aerosol_met_prior_creator:create_parent_object(sub_object)
   -- Luabind can only handle up to 10 arguments per function. As an easy
   -- work around we put the various thresholds into an array.
   local mq = Blitz_double_array_1d(7)
   mq:set(0, self.exp_aod)
   mq:set(1, self.min_types)
   mq:set(2, self.max_types)
   if(self.linear_aod) then
      mq:set(3, 1)
   else
      mq:set(3, 0)
   end
   if(self.relative_humidity_aerosol) then
      mq:set(4, 1)
   else
      mq:set(4, 0)
   end
   mq:set(5, self.max_residual)

   self.aerosol_met_prior = AerosolMetPrior.create(
     self.config.met, self.config:h_merra_aerosol(),
     self.config.pressure,
     self.config.relative_humidity,
     -- In production, have a single covariance for all merra types,
     -- but allow for using functions that have a different covariance
     -- for different types
     self:covariance(0), 
     mq)
   self.sub_object_save = sub_object
   local i, t
   for i, t in ipairs(sub_object) do
      self.aerosol_met_prior:add_aerosol(t.extinction, t.property, t.initial_guess)
   end
   local res = self.aerosol_met_prior:aerosol()
   self.config.number_aerosol = res:number_particle()
   return res
end

function ConfigCommon.aerosol_met_prior_creator:initial_guess()
   local cig = CompositeInitialGuess()
   local flag = ConfigCommon.create_flag(3, Range.all())
   local cig_merra = CompositeInitialGuess()
   cig_merra:add_builder(self.aerosol_met_prior:initial_guess())

   --- Use AOD from merra, but leave the remaining shape parameters at values
   --- fixed at the supplied apriori values (e.g., read from HDF file, or 
   --- whatever).
   local ig_merra = cig_merra:initial_guess()
   for i=1,self.aerosol_met_prior:number_merra_particle() do
      -- Send particle number to ap and cov functions in case they
      -- can use them
      local oval = self:apriori(i-1)

      -- Optionally ignore the merra initial guess value and just use
      -- value from apriori function
      if (self.ignore_merra_aod == nil) then
        oval:set(0, ig_merra((i - 1) * 3))
      end

      local ig = InitialGuessValue()
      ig:apriori_subset(flag, oval)
      ig:apriori_covariance_subset(flag, self:covariance(i-1))
      cig:add_builder(ig)
   end

   --- Now add in the fixed particles.
   local i, t
   for i, t in ipairs(self.sub_object_save) do
      cig:add_builder(t.initial_guess)
   end

   return cig
end

function ConfigCommon.aerosol_met_prior_creator:register_output(ro)
   local h = HdfFile(self.config.met_file)
   local all_aer_names = h:read_string_vector("/Metadata/CompositeAerosolTypes")
   for i, aer_name in ipairs(self.aerosols) do
       all_aer_names:push_back(aer_name)
   end
 
   ro:push_back(AerosolConsolidatedOutput(self.config.aerosol, all_aer_names))
end

------------------------------------------------------------
--- Rayleigh only, we have a nil Aerosol class.
------------------------------------------------------------

ConfigCommon.rayleigh_only = Creator:new()

function ConfigCommon.rayleigh_only:create()
   self.config.number_aerosol = 0
   return nil
end

------------------------------------------------------------
--- Create hydrostatic altitude.
------------------------------------------------------------

ConfigCommon.hydrostatic_altitude = Creator:new()

function ConfigCommon.hydrostatic_altitude:create()
   local res = VectorAltitude()
   local c = self.config
   for i=0,c.l1b:number_spectrometer() - 1 do
      local alts = AltitudeHydrostatic(c.pressure, c.temperature, 
               c.l1b:latitude(i), c.l1b:altitude(i))
      res:push_back(alts)
   end
   return res
end

function ConfigCommon.hydrostatic_altitude:register_output(ro)
   ro:push_back(AltitudeOutput(self.config.altitude:value(0), self.config.pressure))
end

------------------------------------------------------------
--- Create relative humidity. Right now only one class that
--- calculates this, but we have all the infrastructure in 
--- place to allow a class hierarchy if we need this in the
--- future
------------------------------------------------------------

ConfigCommon.calc_relative_humidity = Creator:new()

function ConfigCommon.calc_relative_humidity:create()
   local res
   local c = self.config
   res = RelativeHumidity(c.absorber, c.temperature, c.pressure)
   return res
end

------------------------------------------------------------
--- Map an array to a VectorDouble
------------------------------------------------------------

function ConfigCommon.to_VectorDouble(arr)
   local res = VectorDouble()
   for ind, v in ipairs(arr) do
      res:push_back(v)
   end
   return res
end

------------------------------------------------------------
--- This opens a given absco file, returning a AbscoHdf.
--- This tries to get the absco path name from several places:
---  1. If abscodir environment variable is defined we get it
---     from that.
---  2. If self.absco_local_path is defined, we try to read
---     from that directory. It is ok for this to fail.
---  3. If self.absco_path is defined, we read that.
--- If none of these work, we throw an exception.
------------------------------------------------------------

function ConfigCommon.open_absco(self, fname, table_scale)
   if(os.getenv("abscodir")) then
      return AbscoHdf(os.getenv("abscodir") .. "/" .. fname, table_scale)
   end
   local absco
   if(self.absco_local_path and
      pcall(function() absco = 
               AbscoHdf(self.absco_local_path .. "/" .. fname, 
			table_scale) end)) then
      return absco
   end
   if(self.absco_path) then
      return AbscoHdf(self.absco_path .. "/" .. fname, table_scale)
   end
   error({code=-1})
end

------------------------------------------------------------
--- Handle when table_scale is an array.
------------------------------------------------------------

function ConfigCommon.open_absco_byspecindex(self, fname, sb, table_scale)
   local tsc = ConfigCommon.to_VectorDouble(table_scale)
   if(os.getenv("abscodir")) then
      return AbscoHdf(os.getenv("abscodir") .. "/" .. fname, sb, tsc)
   end
   local absco
   if(self.absco_local_path and
      pcall(function() absco = 
               AbscoHdf(self.absco_local_path .. "/" .. fname, sb, tsc) end)) then
      return absco
   end
   if(self.absco_path) then
      return AbscoHdf(self.absco_path .. "/" .. fname, sb, 
		      tsc)
   end
   error({code=-1})
end

------------------------------------------------------------
--- Common stuff when creating an absorber VMR.
------------------------------------------------------------

CreatorVmr = CreatorApriori:new()

function CreatorVmr:create()
   local res = {}
   if(self.table_scale == nil) then
      res.absco = self.config:open_absco(self.absco, 1.0)
   else
      if(type(self.table_scale) == "number") then
	 res.absco = self.config:open_absco(self.absco, self.table_scale)
      else
	 res.absco = self.config:open_absco_byspecindex(self.absco, 
					self.config:spectral_bound(),
					self.table_scale)
      end
   end
   res.vmr = self:create_vmr()
   return res
end

------------------------------------------------------------
--- Create an absorber fixed level, where we fit for all the
--- VMR levels
------------------------------------------------------------

ConfigCommon.vmr_fixed_level = CreatorVmr:new()

function ConfigCommon.vmr_fixed_level:create_vmr()
   self.vmr = AbsorberVmrFixedLevel(self.config.pressure, self.config.pinp,
                                    self:retrieval_flag(), self:apriori_v(), 
                                    self.name)
   return self.vmr
end

function ConfigCommon.vmr_fixed_level:register_output(ro)
-- Don't write out output if we aren't fitting anything
   local i
   local r
   r = self:retrieval_flag()
   local any_true = false
   for i=1,r:rows() do
      if(r(i - 1)) then
         any_true = true
      end
   end
   if(any_true) then
      ro:push_back(AbsorberVmrFixedLevelOutput.create(self.vmr, 
                                                      self.config.state_vector))
   end
end

------------------------------------------------------------
--- Create an absorber level, where we fit for all the
--- VMR levels
------------------------------------------------------------

ConfigCommon.vmr_level = CreatorVmr:new()

function ConfigCommon.vmr_level:create_vmr()
   self.vmr = AbsorberVmrLevel(self.config.pressure, self:apriori_v(), 
                               self:retrieval_flag(), self.name)
   return self.vmr
end

function ConfigCommon.vmr_level:register_output(ro)
-- Don't write out output if we aren't fitting anything
   local i
   local r
   r = self:retrieval_flag()
   local any_true = false
   for i=1,r:rows() do
      if(r(i - 1)) then
         any_true = true
      end
   end
   if(any_true) then
      ro:push_back(AbsorberVmrLevelOutput.create(self.vmr))
   end
end

------------------------------------------------------------
--- Create an absorber level, where we fit for all the
--- log(VMR) levels
------------------------------------------------------------

ConfigCommon.vmr_log_level = CreatorVmr:new()

function ConfigCommon.vmr_log_level:create_vmr()
   self.vmr = AbsorberVmrLogLevel(self.config.pressure, self:apriori_v(), 
				  self:retrieval_flag(), self.name)
   return self.vmr
end

function ConfigCommon.vmr_log_level:register_output(ro)
-- Don't write out output if we aren't fitting anything
   local i
   local r
   r = self:retrieval_flag()
   local any_true = false
   for i=1,r:rows() do
      if(r(i - 1)) then
         any_true = true
      end
   end
   if(any_true) then
      ro:push_back(AbsorberVmrLogLevelOutput.create(self.vmr))
   end
end

------------------------------------------------------------
--- Create an absorber fixed level, where we hold all the
--- fit for a scale.
------------------------------------------------------------

ConfigCommon.vmr_fixed_level_scaled = CreatorVmr:new()

function ConfigCommon.vmr_fixed_level_scaled:apriori_v()
   local r = Blitz_double_array_1d(1)
   r:set(0, function_or_simple_value(self.scale_apriori, self))
   return r
end

function ConfigCommon.vmr_fixed_level_scaled:covariance_v()
   local r = Blitz_double_array_2d(1, 1)
   r:set(0, 0, function_or_simple_value(self.scale_cov, self))
   return r
end

function ConfigCommon.vmr_fixed_level_scaled:create_vmr()
   self.vmr = AbsorberVmrFixedLevelScaled(self.config.pressure,
                                          self.config.pinp, 
                                          self:apriori(),
                                          self:retrieval_flag()(0),
                                          function_or_simple_value(self.scale_apriori, self), self.name)
   return self.vmr
end

function ConfigCommon.vmr_fixed_level_scaled:register_output(ro)
   ro:push_back(AbsorberVmrFixedLevelScaledOutput.create(self.vmr))
end

------------------------------------------------------------
--- Create an absorber ECMWF, where we fit for a scale.
------------------------------------------------------------

ConfigCommon.vmr_met = CreatorVmr:new()

function ConfigCommon.vmr_met:apriori_v()
   local r = Blitz_double_array_1d(1)
   r:set(0, function_or_simple_value(self.scale_apriori, self))
   return r
end

function ConfigCommon.vmr_met:covariance_v()
   local r = Blitz_double_array_2d(1, 1)
   r:set(0, 0, function_or_simple_value(self.scale_cov, self))
   return r
end

function ConfigCommon.vmr_met:create_vmr()
   self.vmr = AbsorberVmrMet(self.config.met,
                             self.config.pressure,
                             function_or_simple_value(self.scale_apriori, self), 
                             self:retrieval_flag()(0),
                             self.name)
   return self.vmr
end

function ConfigCommon.vmr_met:register_output(ro)
   ro:push_back(AbsorberVmrMetOutput.create(self.vmr))
end

------------------------------------------------------------
--- Create an absorber level, where we fit for a scale.
------------------------------------------------------------

ConfigCommon.vmr_level_scaled = CreatorVmr:new()

function ConfigCommon.vmr_level_scaled:apriori_v()
   local r = Blitz_double_array_1d(1)
   r:set(0, function_or_simple_value(self.scale_apriori, self))
   return r
end

function ConfigCommon.vmr_level_scaled:covariance_v()
   local r = Blitz_double_array_2d(1, 1)
   r:set(0, 0, function_or_simple_value(self.scale_cov, self))
   return r
end

function ConfigCommon.vmr_level_scaled:create_vmr()
   self.vmr = AbsorberVmrLevelScaled(self.config.pressure,
                                     self:vmr_profile(), 
                                     function_or_simple_value(self.scale_apriori, self), 
                                     self:retrieval_flag()(0),
                                     self.name)
   return self.vmr
end

function ConfigCommon.vmr_level_scaled:register_output(ro)
   ro:push_back(AbsorberVmrLevelScaledOutput.create(self.vmr))
end

------------------------------------------------------------
--- Create an absorber fixed level, where we hold all the
--- VMR levels constant
------------------------------------------------------------

ConfigCommon.vmr_fixed_level_constant = ConfigCommon.vmr_fixed_level:new()

function ConfigCommon.vmr_fixed_level_constant:retrieval_flag()
   local res = Blitz_bool_array_1d(self:apriori_v():rows())
   res:set(Range.all(), false)
   return res
end

function ConfigCommon.vmr_fixed_level_constant:initial_guess()
   return CompositeInitialGuess()
end

------------------------------------------------------------
--- Create an absorber level, where we hold all the
--- VMR levels constant
------------------------------------------------------------

ConfigCommon.vmr_level_constant = ConfigCommon.vmr_level:new()

function ConfigCommon.vmr_level_constant:retrieval_flag()
   local res = Blitz_bool_array_1d(self:apriori_v():rows())
   res:set(Range.all(), false)
   return res
end

function ConfigCommon.vmr_level_constant:initial_guess()
   return CompositeInitialGuess()
end

------------------------------------------------------------
--- Create an absorber level, where we hold all the
--- VMR levels constant AND the gas is well mixed
--- meaning that all levels have the same value
------------------------------------------------------------

ConfigCommon.vmr_level_constant_well_mixed = ConfigCommon.vmr_level_constant:new()

function ConfigCommon.vmr_level_constant_well_mixed:apriori_v()
   --- It doesn't really matter the size of the profile since
   --- all levels will have the same value. But just make this
   --- profile consistent with all others.
   const_val = self:apriori(self.name)
   local const_arr = Blitz_double_array_1d(self.config.pressure:max_number_level())
   const_arr:set(Range(), const_val(0))
   return const_arr
end

------------------------------------------------------------
--- Create an absorber fixed level, where we hold all the
--- VMR levels constant AND the gas is well mixed
--- meaning that all levels have the same value
------------------------------------------------------------

ConfigCommon.vmr_fixed_level_constant_well_mixed = ConfigCommon.vmr_fixed_level_constant:new()

function ConfigCommon.vmr_fixed_level_constant_well_mixed:apriori_v()
   --- It doesn't really matter the size of the profile since
   --- all levels will have the same value. But just make this
   --- profile consistent with all others.
   const_val = self:apriori(self.name)
   local const_arr = Blitz_double_array_1d(self.config.pressure:max_number_level())
   const_arr:set(Range(), const_val(0))
   return const_arr
end

------------------------------------------------------------
--- Create absorber.
------------------------------------------------------------

ConfigCommon.absorber_creator = CompositeCreator:new()

function ConfigCommon.absorber_creator:sub_object_key()
   return self.gases
end

function ConfigCommon.absorber_creator:create_parent_object(sub_object)
   local vmr = VectorAbsorberVmr()
   local absco = VectorGasAbsorption()
   local i, t
   for i,t in ipairs(sub_object) do
      vmr:push_back(t.vmr)
      absco:push_back(t.absco)
   end
   if self.number_sub_layers then
      return AbsorberAbsco(vmr, self.config.pressure, self.config.temperature, 
                           self.config.altitude, absco, self.config.constants,
                           self.number_sub_layers)
   else
      return AbsorberAbsco(vmr, self.config.pressure, self.config.temperature, 
                           self.config.altitude, absco, self.config.constants)
   end
end

function ConfigCommon.absorber_creator:register_output(ro)
   CompositeCreator.register_output(self, ro)
   ro:push_back(AbsorberAbscoOutput.create(self.config.absorber, self.config:spectral_bound()))

   -- Create output for GasVmrApriori object it it was created
   -- not really any better places to put this
   if (self.config.ref_co2_ap_obj) then
       ro:push_back(GasVmrAprioriOutput(self.config.ref_co2_ap_obj))
   end
end

------------------------------------------------------------
--- Standard spectral window, from reading the HDF file
------------------------------------------------------------

ConfigCommon.spectral_window_hdf = Creator:new()

function ConfigCommon.spectral_window_hdf:create()
    -- By default read window ranges from static input file, but allow
    -- an option to be specified for supplying a callback to get the
    -- actual ranges.
    local win_ranges
    if (self.window_ranges and self.window_ranges()) then
        win_ranges = self.window_ranges()
    else
        win_ranges = self.config:h():read_double_with_unit_3d("Spectral_Window/microwindow")
    end

    -- If a bad sample mask attribute is supplied for the creator then
    -- use that data when creating the spectral window class
    if (self.bad_sample_mask) then
        local bsamp = self:bad_sample_mask()
        if(bsamp) then
           return SpectralWindowRange(win_ranges, bsamp)
        else
           return SpectralWindowRange(win_ranges)
	end
    else
        return SpectralWindowRange(win_ranges)
    end
end

------------------------------------------------------------
--- Short cut for creating a spectral_bound, including the
--- dispersion needed to convert from sample_index
------------------------------------------------------------

function ConfigCommon:spectral_bound()
   local i
   local d = VectorDispersion()
   for i=1,self.number_pixel:rows() do
      d:push_back(self.dispersion[i])
   end
   self.spec_win:dispersion(d)
   return self.spec_win:spectral_bound()
end

------------------------------------------------------------
--- Base for spectrum sampling creators that does the
--- work of creating an ArrayWithUnit from the high
--- resolution spacing value supplied in the config
--- so that either a single value or a table of values
--- with one per spectrometer.
------------------------------------------------------------

ConfigCommon.spectrum_sampling_base = Creator:new()

function ConfigCommon.spectrum_sampling_base:high_res_spacing()
   local spacing_val = self.high_resolution_spectrum_spacing
   local spacing_awu
   if type(spacing_val) == "table" then
      spacing_arr = Blitz_double_array_1d(#spacing_val)
      for idx,dbl_w_u in ipairs(spacing_val) do
         spacing_arr:set(idx-1, dbl_w_u.value)
      end
      spacing_awu = ArrayWithUnit_1d(spacing_arr, spacing_val[1].units)
   else
      nspec = self.config.number_pixel:rows()
      spacing_arr = Blitz_double_array_1d(nspec)
      spacing_arr:set(Range.all(), spacing_val.value)
      spacing_awu = ArrayWithUnit_1d(spacing_arr, spacing_val.units)
   end
   return spacing_awu
end

------------------------------------------------------------
--- Create uniform spectral sampling
------------------------------------------------------------

ConfigCommon.uniform_spectrum_sampling = ConfigCommon.spectrum_sampling_base:new()

function ConfigCommon.uniform_spectrum_sampling:create()
   return SpectrumSamplingFixedSpacing(self:high_res_spacing())
end

------------------------------------------------------------
--- Nonuniform spectral sampling
---
--- This depends on:
---  self.nonunif_rt_grid_files.o2
---  self.nonunif_rt_grid_files.weak_co2
---  self.nonunif_rt_grid_files.strong_co2
------------------------------------------------------------

ConfigCommon.nonuniform_spectrum_sampling = ConfigCommon.spectrum_sampling_base:new()

function ConfigCommon.nonuniform_spectrum_sampling:create()
   local uspec_samp = SpectrumSamplingFixedSpacing(self:high_res_spacing())
   self.nonunif_rt_grid_files.config = self.config
   return NonuniformSpectrumSampling(self.nonunif_rt_grid_files:o2(), 
                                     self.nonunif_rt_grid_files:weak_co2(), 
                                     self.nonunif_rt_grid_files:strong_co2(),
                                     uspec_samp)
end


------------------------------------------------------------
-- Create Stokes coefficient as constant values (not dependent
-- on state vector).
------------------------------------------------------------

ConfigCommon.stokes_coefficient_constant = Creator:new {}

function ConfigCommon.stokes_coefficient_constant:create()
   return StokesCoefficientConstant(self:value())
end

------------------------------------------------------------
-- Create Stokes coefficient as constant values with a 
-- fit for the fraction of parallel polarization
------------------------------------------------------------

ConfigCommon.stokes_coefficient_fraction = CreatorApriori:new()

function ConfigCommon.stokes_coefficient_fraction:create()
   return StokesCoefficientFraction(self:value(), self:apriori_v(), 
				    self:retrieval_flag())
end

function ConfigCommon.stokes_coefficient_fraction:add_to_statevector(sv)
   sv:add_observer(self.config.stokes_coefficient)
end

function ConfigCommon.stokes_coefficient_fraction:register_output(ro)
   for i=1,self.config.number_pixel:rows() do
      local hdf_band_name = self.config.common.hdf_band_name:value(i-1)
      ro:push_back(StokesCoefficientFractionOutput.create(self.config.stokes_coefficient, i -1, hdf_band_name))
   end
end


------------------------------------------------------------
-- Get stokes coefficients using l1b_stokes_coef. 
------------------------------------------------------------

function ConfigCommon.stokes_coefficient_l1b(self)
   return self.config.l1b:stokes_coef()
end

------------------------------------------------------------
-- Get stokes coefficients using passed in values, as a
-- lua array.
------------------------------------------------------------

function ConfigCommon.stokes_coefficient_value(val)
   return function(self)
	     return ConfigCommon.lua_to_blitz_double_2d(val)
	  end
end

------------------------------------------------------------
--- Just a base Creator class to encapsulate the nadir
--- threshold determination
------------------------------------------------------------

RtCreator = Creator:new()

function RtCreator:pure_nadir()
   local zen = self.config.l1b:zen()

   local nadir_thresh = self.nadir_threshold
   if(nadir_thresh == nil) then
       -- Use old constant value if not defined
       nadir_thresh = 1.0e-6
   end

   -- Do not use pure nadir mode if any of the sounding
   -- zenith angles are above or equal to the threshold
   pure_nadir = true
   for i = 0, zen:rows()-1 do
       if (zen(i) >= nadir_thresh) then
           pure_nadir = false
       end
   end

   return pure_nadir
end

------------------------------------------------------------
--- Standard way of creating RT, made up of LIDORT, LRad
--- and LSI.
------------------------------------------------------------
ConfigCommon.radiative_transfer_lsi = RtCreator:new()

function ConfigCommon.radiative_transfer_lsi:create()
   local sza = self.config.l1b:sza()
   local azm = self.config.l1b:azm()
   local zen = self.config.l1b:zen()
   local do_multiple_scattering_only = true
   local low_stream = self.lsi_constant.low_stream
   local high_stream = self.lsi_constant.high_stream
   local nmom_low = low_stream * 2
   local pure_nadir = self:pure_nadir()

   -- Do not use dedicated twostream solver by default
   -- unless specifically enabled
   local use_twostream
   if(self.lsi_constant.dedicated_twostream == nil) then
       use_twostream = false
   else
       use_twostream = self.lsi_constant.dedicated_twostream
   end

   -- Always use full full quadratures with twostream
   -- unless comparing to LIDORT
   local do_full_quadrature
   if(self.lsi_constant.twostream_full_quadrature == nil) then
       do_full_quadrature = true 
   else
       do_full_quadrature = self.lsi_constant.twostream_full_quadrature
   end

   -- Minimum nmom allowed by LIDORT is 3
   if(nmom_low < 3) then nmom_low = 3 end

   local rt_low
   if(low_stream == 1 and use_twostream) then
      rt_low = TwostreamRt(self.config.atmosphere, self.config.stokes_coefficient,
                           sza, zen, azm, do_full_quadrature)
   else
      rt_low = LidortRt(self.config.atmosphere, self.config.stokes_coefficient,
                        sza, zen, azm, pure_nadir, low_stream, nmom_low, 
                        do_multiple_scattering_only)
   end
   local rt_high = LidortRt(self.config.atmosphere, self.config.stokes_coefficient,
                            sza, zen, azm, pure_nadir, high_stream, LidortRt.maxmoments_input, 
                            do_multiple_scattering_only)
   rt_low = LRadRt.create(rt_low, self.config:spectral_bound(), 
                              sza, zen, azm, pure_nadir, true, false)
   rt_high = LRadRt.create(rt_high, self.config:spectral_bound(), 
                               sza, zen, azm, pure_nadir, true, true)
   rt_high = HresWrapper.create(rt_high)
   return LsiRt.create(rt_low, rt_high, self.config:h(), "LSI")
end

------------------------------------------------------------
--- RT that doesn't have LSI, instead just LRad
---
--- Depends on
---   self.atmosphere
---   self.l1b
---   self.nstream
---   self.spec_win
------------------------------------------------------------

ConfigCommon.radiative_transfer_lrad = RtCreator:new()

function ConfigCommon.radiative_transfer_lrad:create()
   local sza = self.config.l1b:sza()
   local azm = self.config.l1b:azm()
   local zen = self.config.l1b:zen()
   local do_multiple_scattering_only = true
   local nstream = self.nstream
   local nmom = nstream * 2
   local pure_nadir = self:pure_nadir()

   -- Minimum nmom allowed by LIDORT is 3
   if(nmom < 3) then nmom = 3 end

   local rt = LidortRt(self.config.atmosphere, self.config.stokes_coefficient,
                       sza, zen, azm, pure_nadir, nstream, nmom, 
                       do_multiple_scattering_only)

   return LRadRt.create(rt, self.config:spectral_bound(), 
                            sza, zen, azm, pure_nadir, true, false)
end

------------------------------------------------------------
--- RT used for uplooking observations
------------------------------------------------------------

ConfigCommon.chapman_boa_rt = Creator:new()

function ConfigCommon.chapman_boa_rt:create()
   return ChapmanBoaRT.create(self.config.atmosphere, self.config.l1b:sza(), self.config:spectral_bound())
end

------------------------------------------------------------
--- Create state vector.
------------------------------------------------------------
ConfigCommon.state_vector_creator = Creator:new()

function ConfigCommon.state_vector_creator:create()
   return StateVector();
end

function ConfigCommon.state_vector_creator:register_output(ro)
   if (self.config.do_retrieval) then
      ro:push_back(StateVectorOutput(self.config.state_vector))
   end
end

------------------------------------------------------------
--- Create forward model.
------------------------------------------------------------

ConfigCommon.oco_forward_model = CompositeCreator:new()

function ConfigCommon.oco_forward_model:sub_object_key()
   return {"common", "spec_win", "input", "stokes_coefficient", 
	   "instrument", "atmosphere",
           "spec_samp", "spectrum_effect", "rt", "state_vector"}
end

function ConfigCommon.oco_forward_model:sub_initial_guess_key()
   -- For historical reasons, we put the instrument part of the state
   -- vector after the atmosphere, even though we create them in the
   -- opposite order. This is completely arbitrary, but matches what
   -- is the "expected" order in the state vector (i.e., what a scientist
   -- looking at the state vector expects, the L2 code doesn't care at
   -- all)
   return {"input", "atmosphere", "instrument", "spec_win", 
           "spec_samp", "rt", "spectrum_effect", "state_vector",
	   "stokes_coefficient"}
end

function ConfigCommon.oco_forward_model:create_parent_object(sub_object)
   return OcoForwardModel(self.config.instrument,
                          self.config.spec_win, self.config.l1b, 
                          self.config.rt, self.config.spec_samp, 
                          self.config.state_vector,
                          self.config.spectrum_effect)
end

function ConfigCommon.oco_forward_model:register_output(ro)
   -- Store typical ForwardModel output
   CompositeCreator.register_output(self, ro)
   ro:push_back(ForwardModelOutput(self.config.forward_model))
   ro:push_back(OcoForwardModelOutput.create(self.config.forward_model))

   -- Optionally save high resolution spectra
   if(self.config.write_high_res_spectra) then
      hr_spec_out = HighResSpectrumOutput()
      hr_spec_out:add_as_observer(self.config.forward_model)
      hr_spec_out:add_as_observer(self.config.rt)
      ro:push_back(hr_spec_out:as_register_output_base())
   end

   -- Store source filenames in output, probably a better place
   -- to put this, but its probably also fine here
   local dataset_name = VectorString()
   local file_name = VectorString()

   -- Add source data files
   for i,var_ds in ipairs({ { "L1BFile", "spectrum_file"},
                            { "ResampledMetFile", "met_file" },
                            { "CO2PriorFile", "co2_pr_file" },
			    { "StaticInput", "static_file"},
			    { "SolarFile", "static_solar_file"},
			    { "AerosolFile", "static_aerosol_file"},
			    { "MerraFile", "merra_file_name"},
			    { "IMAPFile", "imap_file"},
			    { "EOFFile", "static_eof_file"},
			    { "AtmosphereFile", "atmosphere_file"},
			    { "RunLog", "runlog_file"},
			 }) do
      ds_name = var_ds[1]
      var_name = var_ds[2]
      f_name = self.config[var_name]

      if f_name then
         dataset_name:push_back(ds_name)
         file_name:push_back(f_name)
      else
         self.config.diagnostic_message("Filename for " .. var_name .. " is not defined") 
      end
   end

   -- Add absco filesself.ignore_merra_ig)
   for i,gas_name in ipairs(self.config.fm.atmosphere.absorber.gases) do
      dataset_name:push_back("AbscoFile" .. gas_name)
      file_name:push_back(self.config.fm.atmosphere.absorber[gas_name].absco)
   end
   
   ro:push_back(SourceFilesOutput("Metadata", dataset_name, file_name))
end

------------------------------------------------------------
-- Create the connor solver we normally use.
--
-- Depends on:
--   self.forward_model
--   self.solver_constant
------------------------------------------------------------

function ConfigCommon:connor_solver(config)
   local cost_func = ForwardModelCostFunction(config.forward_model)
   local conv = ConnorConvergence(config.forward_model, 
                                  self.threshold, 
                                  self.max_iteration, 
                                  self.max_divergence, 
                                  self.max_chisq)
   local out = ConnorConvergenceOutput.create(conv)
   config.register_output:push_back(out)
   config.conn_solver = ConnorSolver(cost_func, conv, 
				     self.gamma_initial)
   local iter_log = SolverIterationLog(config.state_vector)
   iter_log:add_as_observer(config.conn_solver)
   out = ConnorSolverOutput(config.conn_solver, config.write_jacobian)
   config.register_output:push_back(out)
end

------------------------------------------------------------
--- Set up error analysis
--- 
--- This depends on:
---   self.solver
---   self.atmosphere
---   self.forward_model
------------------------------------------------------------

function ConfigCommon:create_error_analysis()
   if(self.iter_solver == nil) then
      self.error_analysis = ErrorAnalysis(self.conn_solver, self.atmosphere, 
					  self.forward_model)
   else
      self.error_analysis = ErrorAnalysis(self.stat_method_map,
					  self.atmosphere, 
					  self.forward_model)
   end
   have_co2 = self.absorber:gas_index("CO2") ~= -1
   local out = ErrorAnalysisOutput(self.error_analysis, self:spec_flag(), 
				   have_co2)
   self.register_output:push_back(out)
end

function ConfigCommon:nlls_max_likelihood(config)
   config.stat_method_map = MaxLikelihoodOCO(config.forward_model)
   config.opt_problem = NLLSMaxLikelihood(config.stat_method_map, true)
   config.opt_problem:parameters(config.initial_guess:initial_guess())
end

function ConfigCommon:nlls_max_a_posteriori(config)
   config.stat_method_map = 
      MaxAPosterioriOCO( config.forward_model, 
                         config.initial_guess:apriori(), 
                         config.initial_guess:apriori_covariance() )
   config.opt_problem = NLLSMaxAPosteriori(config.stat_method_map, true)
   config.opt_problem:parameters(config.initial_guess:initial_guess())

   local out = MaxAPosterioriOutput(config.stat_method_map, config.write_jacobian)
   config.register_output:push_back(out)
end


function ConfigCommon:scaled_nlls_max_a_posteriori(config)
   ConfigCommon.nlls_max_a_posteriori(self,config)
   config.opt_problem = NLLSProblemScaled.create(
      config.stat_method_map:param_a_priori_uncertainty(),
      config.opt_problem)
   config.opt_problem:parameters(config.opt_problem:scale_parameters(config.initial_guess:initial_guess()))
end

function ConfigCommon:nlls_solver_gsl_lmsder(config)
   config.iter_solver = 
      NLLSSolverGSLLMSDER.create(self.max_cost_function_calls,
				 self.dx_tol_abs,
				 self.dx_tol_rel,
				 self.g_tol_abs,
				 config.opt_problem,
				 true)
end

function ConfigCommon:iterative_solver(config)
   self:opt_problem(config)
   self:iter_solver(config)
end

-- Replaced with iterative_solver, but leave this here for a bit
--- for reference.
function ConfigCommon:create_problem_and_solver()
--    t = self:create_nlls_max_likelihood()
   self:create_nlls_max_a_posteriori()
--    t = self:create_scaled_nlls_max_a_posteriori()

   --  Retrieval method 1
    self.iterative_solver = NLLSSolverGSLLMSDER.create( self.retrieval_thresholds.max_cost_function_calls,
                                        self.retrieval_thresholds.dx_tol_abs,
                                        self.retrieval_thresholds.dx_tol_rel,
                                        self.retrieval_thresholds.g_tol_abs,
                                        self.opt_problem,
                                        true )
   --  Retrieval method 2
--   t.iter_solver = NLLSSolverGSLLMDER.create( self.retrieval_thresholds.max_cost_function_calls,
--                                       self.retrieval_thresholds.dx_tol_abs,
--                                       self.retrieval_thresholds.dx_tol_rel,
--                                       self.retrieval_thresholds.g_tol_abs,
--                                       t.opt_problem,
--                                       true )

   --  Retrieval method 3
--   local conv = ConnorConvergence( self.forward_model, 
--                                   self.solver_constant.threshold, 
--                                   self.solver_constant.max_iteration, 
--                                   self.solver_constant.max_divergence, 
--                                   self.solver_constant.max_chisq )
--   t.iter_solver = ConnorSolverMAP.create( self.retrieval_thresholds.max_cost_function_calls,
--                                           self.retrieval_thresholds.dx_tol_abs,
--                                           self.retrieval_thresholds.dx_tol_rel,
--                                           self.retrieval_thresholds.g_tol_abs,
--                                           t.opt_problem,
--                                           conv,
--                                           self.solver_constant.gamma_initial )

   --  Retrieval method 4
--   t.iter_solver = CostMinimizerGSL.create( self.retrieval_thresholds.max_cost_function_calls,
--                                     self.retrieval_thresholds.dx_tol_abs,
--                                     self.retrieval_thresholds.dx_tol_rel,
--                                     self.retrieval_thresholds.minimizer_size_tol,
--                                     t.opt_problem,
--                                     t.stat_method_map:param_a_priori_uncertainty()/8.0,
--                                     true )
end

--- Base directory, used to find the solar model
config_common_dir = ConfigCommon.local_dir()
