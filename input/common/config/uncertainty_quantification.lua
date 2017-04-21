-- For table.contains
require "helper_functions"

function init_uq(config)

    config.diagnostic = true

    ------------
    -- Common --
    ------------
    
    function uq_apriori(field)
        return function(self)
            return self.config:l1b_hdf_file():read_double_1d(field .. "/a_priori")
        end
    end

    function uq_apriori_i(field)
        return function(self, i)
            return self.config:l1b_hdf_file():read_double_2d(field .. "/a_priori")(i, Range.all())
        end
    end

    function uq_apriori_sounding_scalar(field)
        return function(self)
            local sid = self.config:l1b_sid_list()
            local frame_idx = sid:frame_number()
            local ap = Blitz_double_array_1d(1)
            ap:set(0, self.config:l1b_hdf_file():read_double_1d(field .. "/a_priori")(frame_idx))
            return ap
        end
    end

    function uq_apriori_sounding(field)
        return function(self)
            local sid = self.config:l1b_sid_list()
            local frame_idx = sid:frame_number()
            return self.config:l1b_hdf_file():read_double_2d(field .. "/a_priori")(Range.all(), frame_idx)
        end
    end

    function uq_apriori_sounding_i(field)
        return function(self, i)
            local sid = self.config:l1b_sid_list()
            local frame_idx = sid:frame_number()
            return self.config:l1b_hdf_file():read_double_3d(field .. "/a_priori")(i, Range.all(), frame_idx)
        end
    end


    function uq_covariance(field)
        return function(self)
            return self.config:l1b_hdf_file():read_double_2d(field .. "/covariance")
        end
    end

    function uq_covariance_i(field)
        return function(self, i)
            return self.config:l1b_hdf_file():read_double_3d(field .. "/covariance")(i, Range.all(), Range.all())
        end
    end

    ---------------
    -- SNR Coefs --
    ---------------

    -- No sounding dimension in UQ files so use a new bad sample mask
    -- function that is aware of this
    function snr_coef_bad_sample_mask_uq(self)
        local l1b_hdf_file = self.config:l1b_hdf_file()
        local snr_coef = l1b_hdf_file:read_double_3d("/InstrumentHeader/snr_coef")

        if (snr_coef:depth() > 2) then
            local sid = self.config:l1b_sid_list()
            local bad_sample_mask = snr_coef(Range.all(), Range.all(), 2)
            return bad_sample_mask
        end
        
        return nil
    end
    function bad_sample_list_bad_sample_mask(self)
        local l1b_hdf_file = self.config:l1b_hdf_file()
	local bad_sample_mask = l1b_hdf_file:read_double_2d("/InstrumentHeader/bad_sample_list")
	return bad_sample_mask
    end
    function l1b_bad_sample_mask(self)
       local l1b_hdf_file = self.config:l1b_hdf_file()
       if(l1b_hdf_file:has_object("/InstrumentHeader/bad_sample_list")) then
	  return self.config.bad_sample_list_bad_sample_mask(self)
       else
	  return self.config.snr_coef_bad_sample_mask(self)
       end
    end   

    config.fm.spec_win.bad_sample_mask = l1b_bad_sample_mask_uq

    -----------
    -- Noise --
    -----------

    -- Redefined noise function due to change in dataset sizes
    local uq_noise = Creator:new()

    function uq_noise:create(hdf_file, sid)
        local hv = self.config:l1b_hdf_file()
        local sid = self.config:l1b_sid_list()
        local nspec = self.config.spec_win:number_spectrometer()
        local max_ms = Blitz_double_array_1d(nspec)
        for i, v in ipairs(self.max_ms) do
            max_ms:set(i-1, v)
        end
        return UqNoiseModel(hv, sid, max_ms)
    end

    config.fm.input.l1b.noise.creator = uq_noise

    ---------
    -- L1B --
    ---------

    -- Redefine the sounding id list function to use a class that
    -- can handle the UQ Sounding Ids
    function config:l1b_sid_list()
        if (not self.sid_string or self.sid_string == "") then
            error("Sounding id string not defined")
        end
        if(not self.l1b_sid_list_v) then
           local hv = self:l1b_hdf_file()
           if(hv) then
              self.l1b_sid_list_v = UqSoundingId(hv, self.sid_string)
              self.sid = self.l1b_sid_list_v
           end
        end
        return self.l1b_sid_list_v
    end

    -- Use custom L1B class
    local level1b_uq = CreatorL1b:new()

    function level1b_uq:create_parent_object()
       local hv = self.config:l1b_hdf_file()
       local sid = self.config:l1b_sid_list()
       local l1b_uq = Level1bUq(hv, sid)
       l1b_uq.noise_model = self.config.noise
       return l1b_uq
    end
 
    config.fm.input.l1b.creator = level1b_uq

    -----------
    -- ECMWF --
    -----------

    uq_met = Creator:new()

    function uq_met:create()
        local met = UqEcmwf(self.config.spectrum_file)
        self.config.input_file_description = self.config.input_file_description .. "ECMWF input file:    " .. self.config.spectrum_file .. "\n"
        return met
    end

    function uq_met:register_output(ro)
        if (self.config.met) then
            ro:push_back(MetPassThroughOutput(self.config.met))
        end
    end

    config.fm.input.met.creator = uq_met

    ---------
    -- ILS --
    ---------

    -- Use custom ILS table creator due to change in dataset shape
    local ils_table_uq = Creator:new()

    function ils_table_uq:create()
        -- Load reference to L1B HDF file
        local l1b_hdf_file = self.config:l1b_hdf_file()
 
        local sid = self.config:l1b_sid_list()
 
        -- Use a simple table ILS since the ILS supplies
        -- a full set of delta_lambda/repsonses per pixel
        interpolate = false
 
        local idx
        local wavenumber, delta_lambda, response_name
        local res = {}
        for idx = 0, self.config.number_pixel:rows() - 1 do
            local hdf_band_name = self.config.common.hdf_band_name:value(idx)
            local desc_band_name = self.config.common.desc_band_name:value(idx)
  
            -- Delta lambda
            delta_lambda = l1b_hdf_file:read_double_3d("/InstrumentHeader/ils_delta_lambda")
            delta_lambda = delta_lambda(idx, Range.all(), Range.all())
  
            -- Load ILS response values
            response = l1b_hdf_file:read_double_3d("/InstrumentHeader/ils_relative_response")
            response = response(idx, Range.all(), Range.all())
  
            -- Calculate center wavenumber from dispersion, there should be number pixel
            -- of these per spectrometer
            wavenumber = self.config.dispersion[idx+1]:pixel_grid():data()
  
            res[idx+1] = IlsTableLog(wavenumber, delta_lambda, response, desc_band_name, hdf_band_name, interpolate)
        end
        return res
    end

    config.fm.instrument.ils_func.creator = ils_table_uq

    ----------------
    -- Dispersion --
    ----------------
    
    -- A priori and covariance from UQ file only supply values for the dispersion offset term
   
    local orig_disp_ap = config.fm.instrument.dispersion.apriori
    function uq_dispersion_ap(self, i)
        local ap = orig_disp_ap(self, i)

        local l1b_hdf_file = self.config:l1b_hdf_file()
        local sid = self.config:l1b_sid_list()
        local frame_idx = sid:frame_number()
        local offset_mean = l1b_hdf_file:read_double_2d("/Instrument/Dispersion/a_priori")(i, frame_idx)

        ap.value:set(0, offset_mean)

        return ap
    end

    local orig_disp_cov = config.fm.instrument.dispersion.covariance
    function uq_dispersion_cov(self, i)
        local cov = orig_disp_cov(self, i)

        local l1b_hdf_file = self.config:l1b_hdf_file()
        local offset_cov = l1b_hdf_file:read_double_2d("/Instrument/Dispersion/covariance")(i, i)

        cov:set(0, 0, offset_cov)

        return cov
    end
    
    config.fm.instrument.dispersion.apriori = uq_dispersion_ap
    config.fm.instrument.dispersion.covariance = uq_dispersion_cov

    ------------
    -- Ground --
    ------------

    -- UQ code only uses lambertian or coxmunk, so modify this since default code now uses BRDF
    function config:uq_ground_type_name()
       if(self:land_or_water() == "land") then
          return "lambertian"
       else
          return "coxmunk"
       end
    end

    config.ground_type_name = config.uq_ground_type_name

    config.fm.atmosphere.ground.lambertian.apriori = uq_apriori_sounding_i("Ground/Albedo")
    config.fm.atmosphere.ground.lambertian.covariance = uq_covariance_i("Ground/Albedo")

    config.fm.atmosphere.ground.coxmunk.apriori = uq_apriori_sounding_scalar("Ground/Windspeed")
    config.fm.atmosphere.ground.coxmunk.covariance = uq_covariance("Ground/Windspeed")

    config.fm.atmosphere.ground.coxmunk_lambertian.apriori = uq_apriori_sounding_i("Ground/Albedo")
    config.fm.atmosphere.ground.coxmunk_lambertian.covariance = uq_covariance_i("Ground/Albedo")

    ---------
    -- Gas --
    ---------

    config.fm.atmosphere.absorber.CO2.apriori = uq_apriori_sounding("Gas/CO2")
    config.fm.atmosphere.absorber.CO2.covariance = uq_covariance("Gas/CO2")

    function uq_h2o_ap(self)
        local l1b_hdf_file = self.config:l1b_hdf_file()
        return l1b_hdf_file:read_double_1d("Gas/H2O_Scaling_factor/a_priori")(0)
    end


    function uq_h2o_cov(self)
        local l1b_hdf_file = self.config:l1b_hdf_file()
        return l1b_hdf_file:read_double_2d("Gas/H2O_Scaling_factor/covariance")(0, 0)
    end

    config.fm.atmosphere.absorber.H2O.scale_apriori = uq_h2o_ap
    config.fm.atmosphere.absorber.H2O.scale_cov = uq_h2o_cov

    -----------------
    -- Temperature --
    -----------------

    config.fm.atmosphere.temperature.apriori = uq_apriori_sounding_scalar("Temperature/Offset")
    config.fm.atmosphere.temperature.covariance = uq_covariance("Temperature/Offset")
    
    ----------------------
    -- Surface Pressure --
    ----------------------

    config.fm.atmosphere.pressure.apriori = uq_apriori_sounding_scalar("Surface_Pressure")
    config.fm.atmosphere.pressure.covariance = uq_covariance("Surface_Pressure")

    -------------
    -- Aerosol --
    -------------
    
    -- Set Ice and Water apriori and covariances
    config.fm.atmosphere.aerosol.Water.apriori = uq_apriori_sounding("Aerosol/Water/Gaussian/Log")
    config.fm.atmosphere.aerosol.Water.covariance = uq_covariance("Aerosol/Water/Gaussian/Log")

    config.fm.atmosphere.aerosol.Ice.apriori = uq_apriori_sounding("Aerosol/Ice/Gaussian/Log")
    config.fm.atmosphere.aerosol.Ice.covariance = uq_covariance("Aerosol/Ice/Gaussian/Log")

    -- Set Merra particle apriori and covariances by modifying the creator to not use Merra
    -- initial guess value
    config.fm.atmosphere.aerosol.ignore_merra_aod = true
    config.fm.atmosphere.aerosol.apriori = uq_apriori_sounding_i("Aerosol/Merra/Gaussian/Log")
    config.fm.atmosphere.aerosol.covariance = uq_covariance_i("Aerosol/Merra/Gaussian/Log")

    ------------------
    -- Fluorescence -- 
    ------------------

    -- Reconfigure portions of the fluoresence creator to not rely on getting data from the L1B radiance routine
    uq_fluorescence = ConfigCommon.fluorescence_effect:new()

    function uq_fluorescence:create()
       local ground_type = self.config:ground_type_name()

       -- Only use this for coxmunk runs
       if(ground_type == "lambertian") then
            local flag = self:retrieval_flag()
 
            local rad_unit = Unit("Ph sec^{-1} m^{-2} sr^{-1} um^{-1}")
            self.fluorescence_effect = FluorescenceEffect(self:apriori(), flag,
                               self.config.atmosphere,
                               self.config.stokes_coefficient,
                               self.config.l1b:zen_with_unit(0),
                               0, self:reference_point(), rad_unit)
            local res = {}
            res[1] =  self.fluorescence_effect
            return res
        else
            return nil
        end
    end

    config.fm.spectrum_effect.fluorescence.creator = uq_fluorescence
    config.fm.spectrum_effect.fluorescence.apriori = uq_apriori_sounding("Fluorescence")
    config.fm.spectrum_effect.fluorescence.covariance = uq_covariance("Fluorescence")

    -- Remove EOF from state vector
    config.fm.instrument.instrument_correction.ic_nadir = {}
    config.fm.instrument.instrument_correction.ic_glint = {}
    config.fm.instrument.instrument_correction.ic_target = {}

end
