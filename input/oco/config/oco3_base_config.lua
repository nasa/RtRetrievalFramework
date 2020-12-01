------------------------------------------------------------
--- OCO-2 is almost the same as OCO-2. Right now, the
--- differences are a different EOF file and a different spectral
--- window

require "oco_base_config"
oco3_base_config_dir = ConfigCommon.local_dir()

Oco3BaseConfig = OcoBaseConfig:new {
   static_file = oco3_base_config_dir .. "/../input/l2_oco3_static_input.h5",

   -- Selected based on L1B2 land or water flag
   static_eof_file_land = oco3_base_config_dir .. "/../input/l2_oco3_eof_ioc_B10208_QTS_land_sorted_L2.h5",
   static_eof_file_ocean = oco3_base_config_dir .. "/../input/l2_oco3_eof_ioc_B10208_QTS_oceanG_alt3_falt1_sorted_L2.h5",

   -- Will get filled in when the above is picked, used for output
   static_eof_file = ""
}

-- Override behavior controlling loading of EOF file, for OCO-3 pick a different file
-- for land than for ocean using the land/water flag
function Oco3BaseConfig:h_eof()
   -- Use static_eof_file if found, otherwise use the same static input
   -- file that we use for everything else (self:h()).
   if(not self.h_eof_v) then
      if(self:land_or_water() == "land") then
          self.static_eof_file = self.static_eof_file_land
      else
          self.static_eof_file = self.static_eof_file_ocean
      end

      self.h_eof_v = HdfFile(self.static_eof_file)
      self.input_file_description = self.input_file_description .. "EOF input file:    " .. self.static_eof_file .. "\n"
   end

   return self.h_eof_v
end
