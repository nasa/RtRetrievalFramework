-- For table.contains
require "helper_functions"

function init_rrv(config)
    -- Controls how big the pressure grid will be
    if (not config.profile_step_size) then
        config.profile_step_size = 106
    end

    -- Read RRV specific data from a RRV file
    config.rrv_file = os.getenv("rrv_file")

    -- Function for reading and storing RRV file
    function r(self)
       if(self.config.rrv_file and not self.r_v) then 
          self.r_v = HdfFile(self.config.rrv_file) 
       end
       return self.r_v
    end

    -- Sample the RRV file's data at the profile step size
    function sample_levels(self, profile)
       local step_size = self.config.profile_step_size
       local res = Blitz_double_array_1d(math.floor(profile:rows() / step_size))
       local out_idx = 0

       -- Start at an index that will allow ending on the final index
       local start_idx = profile:rows() % step_size + self.config.profile_step_size - 1

       for in_idx = start_idx,profile:rows(),step_size do
          res:set(out_idx, profile(in_idx))
          out_idx = out_idx + 1
       end
       return res
    end

    --------------
    -- Pressure --
    --------------

    function read_rrv_pressure(self)
        -- hPa to Pa
        local pip = r(self):read_double_1d("/radiosonde/sonde_pressure") * 100
        return pip
    end

    function rrv_pressure_sigma_a(self)
       if (r(self):has_object("/radiosonde/sonde_pressure")) then
           local press_levels = read_rrv_pressure(self)
           local res = Blitz_double_array_1d(math.floor(press_levels:rows() / self.config.profile_step_size))
           res:set(Range.all(), 0.0)
           return res
       else
           return self.config:h():read_double_1d("Pressure/Pressure_sigma_a")
       end 
    end

    function rrv_pressure_sigma_b(self)
       if (r(self):has_object("/radiosonde/sonde_pressure")) then
           local press_levels = read_rrv_pressure(self)
           local sig_b = sample_levels(self, press_levels)
           sig_b = sig_b / sig_b(sig_b:rows()-1)
           self.config:diagnostic_message("Sigma b = ", sig_b)
           return sig_b
       else
           return self.config:h():read_double_1d("Pressure/Pressure_sigma_b")
       end
    end

    function read_rrv_psurf(self)
       -- hPa to Pa
       local psurf = Blitz_double_array_1d(1)
       psurf:set(0, r(self):read_double_scalar("/meteorological/met_pressure") * 100)
       self.config:diagnostic_message("Surface pressure", psurf)
       return psurf
    end

    function psurf_cov_no_ret(self)
       -- Make covariance real small so is not really retrieved,
       -- work around auto derivative mismatch if psurf not included in state vector
       local cov = ConfigCommon.hdf_covariance("Surface_Pressure")(self) 
       cov:set(0, 0, 1e-20)
       return cov
    end

    config.fm.atmosphere.pressure = {
       apriori = read_rrv_psurf,
       covariance = psurf_cov_no_ret,
       a = rrv_pressure_sigma_a,
       b = rrv_pressure_sigma_b,
       creator = ConfigCommon.pressure_sigma,
    }

    -----------------
    -- Temperature --
    -----------------

    function read_rrv_temperature(self)
       -- Read in Kelvin
       local temp = sample_levels(self, r(self):read_double_1d("/radiosonde/sonde_temperature"))
       self.config:diagnostic_message("Temperature ", temp)
       return temp
    end

    ------------------------------------------------------------
    --- Use ECMWF or RRV temperature data depending if sounde
    --- data is present
    ------------------------------------------------------------

    temperature_rrv_or_met = CreatorApriori:new {}

    function temperature_rrv_or_met:create()
        if (r(self):has_object("/radiosonde/sonde_temperature")) then
           return TemperatureLevelOffset(self.config.pressure, self:temperature_levels(),
                                         self:apriori()(0), self:retrieval_flag()(0))
        else
           return TemperatureMet(self.config.met, self.config.pressure,
                                   self:apriori()(0), self:retrieval_flag()(0))
        end
    end

    function ConfigCommon.temperature_met:register_output(ro)
        if (r(self):has_object("/radiosonde/sonde_temperature")) then
            ro:push_back(TemperatureLevelOffsetOutput.create(self.config.temperature))
        else
            ro:push_back(TemperatureMetOutput.create(self.config.temperature))
        end
    end


    config.fm.atmosphere.temperature = {
       temperature_levels = read_rrv_temperature,
       apriori = ConfigCommon.hdf_apriori("Temperature/Offset"),
       covariance = ConfigCommon.hdf_covariance("Temperature/Offset"),
       creator = temperature_rrv_or_met,
       retrieved = false,
    }

    ---------
    -- H2O --
    ---------

    function read_rrv_h2o(self)
       local h2o = sample_levels(self, r(self):read_double_1d("/radiosonde/sonde_vmr"))
       self.config:diagnostic_message("H2O VMR: ", h2o)
       return h2o
    end

    vmr_rrv_or_met = CreatorVmr:new {}

    function vmr_rrv_or_met:apriori_v()
       local r = Blitz_double_array_1d(1)
       r:set(0, self.scale_apriori)
       return r
    end

    function vmr_rrv_or_met:covariance_v()
       local r = Blitz_double_array_2d(1, 1)
       r:set(0, 0, self.scale_cov)
       return r
    end

    function vmr_rrv_or_met:create_vmr()
        if(r(self):has_object("/radiosonde/sonde_vmr")) then
            self.vmr = AbsorberVmrLevelScaled(self.config.pressure,
                                             self:vmr_profile(), 
                                             self.scale_apriori, 
                                             self:retrieval_flag()(0),
                                             self.name)
        else
            self.vmr = AbsorberVmrMet(self.config.met,
                                        self.config.pressure,
                                        self.scale_apriori, 
                                        self:retrieval_flag()(0),
                                        self.name)
        end
        return self.vmr
    end

    function vmr_rrv_or_met:register_output(ro)
        if(r(self):has_object("/radiosonde/sonde_vmr")) then
            ro:push_back(AbsorberVmrLevelScaledOutput.create(self.vmr))
        else
            ro:push_back(AbsorberVmrMetOutput.create(self.vmr))
        end
    end

    config.fm.atmosphere.absorber.H2O.creator = vmr_rrv_or_met
    config.fm.atmosphere.absorber.H2O.vmr_profile = read_rrv_h2o
    config.fm.atmosphere.absorber.H2O.retrieved = false

    ---------
    -- CO2 --
    ---------

    -- Make CO2 present but not retrieved as it needs to be present for
    -- averaging kernel stuff

    -- This covariance is only used if the profile is retrieved
    function co2_cov(self)
       local sig_a = rrv_pressure_sigma_a(self)

       local res = Blitz_double_array_2d(sig_a:rows(),sig_a:rows())
       res:set(Range.all(), Range.all(), 0.0)

       for idx = 0, res:rows()-1 do
          res:set(idx,idx, 3e-09)
       end

       return res
    end
    config.fm.atmosphere.absorber.CO2.creator = ConfigCommon.vmr_level_scaled
    config.fm.atmosphere.absorber.CO2.vmr_profile = config.fm.atmosphere.absorber.CO2.apriori
    config.fm.atmosphere.absorber.CO2.covariance = co2_cov
    config.fm.atmosphere.absorber.CO2.scale_apriori = 1.0
    config.fm.atmosphere.absorber.CO2.scale_cov = 1e-20
    config.fm.atmosphere.absorber.CO2.retrieved = true 

    ------------
    -- Ground --
    ------------

    -- Redefine so that we always use land, since we know our RRV data is over land
    always_lambertian = DispatchCreator:new()
    function always_lambertian:get_creator()
       return ConfigCommon.ground_lambertian
    end

    config.fm.atmosphere.ground.creator = always_lambertian

    function rrv_lambertian_ap(self, idx)
       local sid = self.config:l1b_sid_list()
       local frame_num = sid:frame_number()
       local sounding_num = sid:sounding_number()

       local refl = r(self):read_double_1d("/surface/ReflCtrWavelength")(idx)
       local offnadir_ratio = r(self):read_double_3d("/sounding/offnadir_nadir_ratio")(frame_num, sounding_num, idx)

       -- insitu file has footprint asd ratios shaped frame_id * sounding_pos
       local fp_asd_ratio = r(self):read_double_3d("/sounding/fp_asd_ratio")(frame_num, sounding_num, idx)

       local ap = Blitz_double_array_1d(2)
       ap:set(0, refl * fp_asd_ratio * offnadir_ratio)
       ap:set(1, 0.0)

       self.config:diagnostic_message(string.format("Albedo (%d): %0.2f * %0.2f * %0.2f = %0.2f", idx, refl, fp_asd_ratio, offnadir_ratio, ap(0)))

       return ap
    end

    config.fm.atmosphere.ground.lambertian.apriori = rrv_lambertian_ap

    orig_lambertian_cov = config.fm.atmosphere.ground.lambertian.covariance
    function rrv_lambertian_cov(self, i)
        local apcov = orig_lambertian_cov(self, i)
        apcov:set(Range.all(), Range.all(), 1e-20)
        return apcov
    end
    config.fm.atmosphere.ground.lambertian.covariance = rrv_lambertian_cov

    rrv_lamberitan_retrieval = ConfigCommon.lambertian_retrieval:new {}

    function rrv_lamberitan_retrieval:create()
       local num_coeff = self:apriori_v(0):rows()
       local num_spec = self.config.number_pixel:rows()

       local ap = Blitz_double_array_2d(num_spec, num_coeff)
       local flag = Blitz_bool_array_2d(num_spec, num_coeff)

       for i = 1, num_spec do
           ap:set(i-1, Range.all(), self:apriori_v(i - 1))
           flag:set(i-1, Range.all(), self:retrieval_flag(i))
       end

       -- Use center wavelength from RRV file
       local center_wl = r(self):read_double_1d("/surface/BandCtrWavelength")
       local band_ref = ArrayWithUnit_1d(center_wl, "micron")

       local lambertian = GroundLambertian(ap, flag, band_ref,
                                           self.config.common.desc_band_name)
       return lambertian
    end

    config.fm.atmosphere.ground.lambertian.creator = rrv_lamberitan_retrieval 
    config.fm.atmosphere.ground.lambertian.retrieve_bands = { false, false, false }

    ----------------
    -- Instrument --
    ----------------

    -- Enable radiance scaling
    config.fm.instrument.instrument_correction.radiance_scaling = {
        apriori = ConfigCommon.hdf_apriori_i("Instrument/RadianceScaling/RRV"),
        covariance = ConfigCommon.hdf_covariance_i("Instrument/RadianceScaling/RRV"),
        creator = ConfigCommon.radiance_scaling_sv_fit,
        retrieve_bands = { true, true, true },
    }

    if (config.fm.instrument.instrument_correction.ic_h_gain) then
        table.insert(config.fm.instrument.instrument_correction.ic_h_gain, 'radiance_scaling')
        table.insert(config.fm.instrument.instrument_correction.ic_m_gain, 'radiance_scaling')
    else
        table.insert(config.fm.instrument.instrument_correction.ic_glint, 'radiance_scaling')
        table.insert(config.fm.instrument.instrument_correction.ic_nadir, 'radiance_scaling')
        table.insert(config.fm.instrument.instrument_correction.ic_target, 'radiance_scaling')
    end

    -- Disable zero offset correction retrieval
    if (config.fm.instrument.instrument_correction.zero_offset_waveform) then
        config.fm.instrument.instrument_correction.zero_offset_waveform.retrieve_bands = { false, false, false }
    end

    -- Disable dispersion retrieval
    config.fm.instrument.dispersion.retrieved = false

    -----------
    -- Other --
    -----------

    -- Disable fluorescence
    if (config.fm.spectrum_effect.speceff_h_gain) then
        -- Remove for GOSAT configs
        if (table.contains(config.fm.spectrum_effect.speceff_h_gain, "fluorescence")) then
            remove_at = table.index(config.fm.spectrum_effect.speceff_h_gain, "fluorescence")
            table.remove(config.fm.spectrum_effect.speceff_h_gain, remove_at)
        end
        if (table.contains(config.fm.spectrum_effect.speceff_m_gain, "fluorescence")) then
            remove_at = table.index(config.fm.spectrum_effect.speceff_m_gain, "fluorescence")
            table.remove(config.fm.spectrum_effect.speceff_m_gain, remove_at)
        end
    else
        if (table.contains(config.fm.spectrum_effect.speceff, "fluorescence")) then
            remove_at = table.index(config.fm.spectrum_effect.speceff, "fluorescence")
            table.remove(config.fm.spectrum_effect.speceff, remove_at)
        end
    end
    config.fm.spectrum_effect.fluorescence.retrieved = false

end
