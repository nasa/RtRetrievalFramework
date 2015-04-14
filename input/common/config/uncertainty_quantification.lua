-- For table.contains
require "helper_functions"

function init_uq(config)

    config.diagnostic = true

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


    -- Use UQ specific reader for ECMWF file
    function config:ecmwf()
        if(self.spectrum_file and not self.ecmwf_v) then 
            self.ecmwf_v = UqEcmwf(self.spectrum_file)
            self.input_file_description = self.input_file_description .. "ECMWF input file:    " .. self.spectrum_file .. "\n"
        end
        return self.ecmwf_v
    end

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

    config.fm.spec_win.bad_sample_mask = snr_coef_bad_sample_mask_uq

    -- Redefined noise function due to change in dataset sizes
    uq_noise = Creator:new()

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

    config.fm.l1b.noise.creator = uq_noise

    -- Use custom L1B class
    level1b_uq = CreatorL1b:new()

    function level1b_uq:create_parent_object()
       local hv = self.config:l1b_hdf_file()
       local sid = self.config:l1b_sid_list()
       local l1b_uq = Level1bUq(hv, sid)
       l1b_uq.noise_model = self.config.noise
       return l1b_uq
    end
 
    config.fm.l1b.creator = level1b_uq

    -- Use custom ILS table creator due to change in dataset shape
    ils_table_uq = Creator:new()

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

    -- Use apriori value from the UQ file for albedo instead of the usual method of deriving it from the radiance
    function albedo_apriori_uq(self, i)
        local l1b_hdf_file = self.config:l1b_hdf_file()
        return l1b_hdf_file:read_double_2d("/Ground/Albedo/a_priori")(i, Range.all())
    end
    config.fm.atmosphere.ground.lambertian.apriori = albedo_apriori_uq

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
    config.fm.spectrum_effect.fluorescence.covariance = ConfigCommon.hdf_covariance("Fluorescence")

    -- Remove EOF from state vector
    config.fm.instrument.instrument_correction.ic_nadir = {}
    config.fm.instrument.instrument_correction.ic_glint = {}
    config.fm.instrument.instrument_correction.ic_target = {}

end
