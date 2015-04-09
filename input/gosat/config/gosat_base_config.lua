-- This sets up the "standard" run we use in the unit tests, you can then
-- use this and override any feature you want to be done differently
--
-- Pull anything from CommonConfig through the CommonConfig name and
-- anything in AcosConfig from that name instead of all through AcosConfig
-- to better differentiate where routines are being loaded from.

require "acos_config"
gosat_base_config_dir = ConfigCommon.local_dir()

GosatBaseConfig = AcosConfig:new {
------------------------------------------------------------
--- Various constants used to describe input data.
------------------------------------------------------------

   -- As a convenience, we pass in a few fields by environment
   -- variables.
   sid_string = os.getenv("sounding_id"),
   spectrum_file = os.getenv("spectrum_file"),
   ecmwf_file = os.getenv("ecmwf_file"),
   static_file = gosat_base_config_dir .. "/../input/l2_gosat_static_input.h5",
   static_solar_file = config_common_dir .. "/../input/l2_solar_model.h5",
   static_aerosol_file = config_common_dir .. "/../input/l2_aerosol_combined.h5",
   merra_dir = "/groups/algorithm/l2_fp/merra_composite",

------------------------------------------------------------
-- Set this true to get diagnostic messages to help debug
-- problems with Lua
------------------------------------------------------------

   diagnostic = false,

------------------------------------------------------------
-- Paths for Absco data. We first look in local path, and if
-- we can't read the file try the full path. This allows us
-- to use a local disk on the cluster.
------------------------------------------------------------

   absco_local_path = "/state/partition1/groups/algorithm/l2_fp/absco",
   absco_path = "/groups/algorithm/l2_fp/absco",

------------------------------------------------------------
-- Connor solver
------------------------------------------------------------

   solver ={ max_iteration=7, max_divergence=2,
	     max_chisq=1.4, threshold=0.1, gamma_initial=5.0,
	     create = ConfigCommon.connor_solver,
	   },

------------------------------------------------------------
-- Iterative solver
--
-- Thresholds and limits for the retrieval process.  
-- Some thresholds are for testing the cost function.  For
-- example g_tol_abs is a tolerance for checking the 
-- gradient of the cost function.  On the other hand, some
-- are used for testing the minimizer (solver).  For
-- example, minimizer_size_tol is used to check the solver.
-- Therefore, we cal this group of constants the retrieval
-- thresholds and not the solver or the problem thresholds.
------------------------------------------------------------

   -- solver = 
   not_used_solver = 
      { max_cost_function_calls=20,
	dx_tol_abs=1e-5, dx_tol_rel=1e-5, 
	g_tol_abs=1e-5,
	minimizer_size_tol=1e-5,
	opt_problem = ConfigCommon.nlls_max_a_posteriori,
	iter_solver = ConfigCommon.nlls_solver_gsl_lmsder,
	create = ConfigCommon.iterative_solver,
     },

------------------------------------------------------------
-- If true then launch solver, otherwise just do a 
-- forward model calculation, with jacobians if
-- write_jacobians is true
------------------------------------------------------------

   do_retrieval = true,

------------------------------------------------------------
-- True if we want to write the jacobians out
------------------------------------------------------------

   write_jacobian = false,

------------------------------------------------------------
-- True if we want to write high resolution spectra out
------------------------------------------------------------

   write_high_res_spectra = false,

------------------------------------------------------------
-- True if we want to generate output for every iteration
------------------------------------------------------------

   iteration_output = false,

------------------------------------------------------------
--- Log level
------------------------------------------------------------

   log_level = LogImp.INFO,

-----------------------------------------------------------
--- Default creators for everything. We default to getting stuff
--- from the HDF file, but this can be changed if desired.
------------------------------------------------------------

   fm = {
      creator = ConfigCommon.oco_forward_model,
      common = {
         desc_band_name = ConfigCommon.hdf_read_string_vector("Common/desc_band_name"),
         hdf_band_name = ConfigCommon.hdf_read_string_vector("Common/hdf_band_name"),
         band_reference = ConfigCommon.hdf_read_double_with_unit_1d("Common/band_reference_point"),
         creator = ConfigCommon.table_function_eval,
      },
      spec_win = {
         creator = ConfigCommon.spectral_window_hdf,
      },
      l1b = {
         creator = AcosConfig.level1b_hdf,
         noise = {
            creator = AcosConfig.gosat_noise_l1b,
         },
      },
      stokes_coefficient = {
	 creator = ConfigCommon.stokes_coefficient_constant,
	 value = ConfigCommon.stokes_coefficient_l1b,
      },
      instrument = {
         creator = ConfigCommon.ils_instrument,
         ils_half_width = { DoubleWithUnit(20, "cm^-1"),
                            DoubleWithUnit(20, "cm^-1"),
                            DoubleWithUnit(20, "cm^-1") },
         dispersion = {
            creator = ConfigCommon.dispersion_polynomial_fitted,
            aband_solar_line_location = DoubleWithUnit(12985.16325e0, "cm^-1"),
            aband_solar_line_width = DoubleWithUnit(0.2*0.2, "cm^-1"),
            aband_search_width = DoubleWithUnit(2.83, "cm^-1"),
            aband_ils_offset = DoubleWithUnit(0.203e0, "cm^-1"),
            offset_scaling = AcosConfig.gosat_offset_scaling("Instrument/Dispersion/offset_scaling"),
            apriori = ConfigCommon.l1b_spectral_coefficient,
            covariance = ConfigCommon.hdf_covariance_i("Instrument/Dispersion"),
            number_pixel = ConfigCommon.hdf_read_int_1d("Instrument/Dispersion/number_pixel"),
            is_one_based = true,
         },
         ils_func = {
            creator = ConfigCommon.ils_table,
         },
         instrument_correction = {
            creator = AcosConfig.instrument_correction_list_h_and_m,
            ic_h_gain = { "zero_offset_waveform", "eof_h_gain_1" },
            ic_m_gain = { "zero_offset_waveform", "eof_m_gain_1" },
            zero_offset_waveform = {
               apriori = ConfigCommon.hdf_apriori_i("Instrument/ZeroLevelOffset"),
               covariance = ConfigCommon.hdf_covariance_i("Instrument/ZeroLevelOffset"),
               creator = ConfigCommon.zero_offset_waveform,
               retrieve_bands = { true, false, false },
            },
            eof_h_gain_1 = {
	       hdf_group = "Instrument/EmpiricalOrthogonalFunction/H_gain",
               apriori = ConfigCommon.hdf_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               covariance = ConfigCommon.hdf_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               order = 1,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = { true, true, true },
            },
            eof_h_gain_2 = {
	       hdf_group = "Instrument/EmpiricalOrthogonalFunction/H_gain",
               apriori = ConfigCommon.hdf_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               covariance = ConfigCommon.hdf_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               order = 2,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = { true, true, true },
            },
            eof_m_gain_1 = {
	       hdf_group = "Instrument/EmpiricalOrthogonalFunction/M_gain",
               apriori = ConfigCommon.hdf_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               covariance = ConfigCommon.hdf_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               order = 1,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = { true, true, true },
            },
            eof_m_gain_2 = {
	       hdf_group = "Instrument/EmpiricalOrthogonalFunction/M_gain",
               apriori = ConfigCommon.hdf_apriori_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               covariance = ConfigCommon.hdf_covariance_i_j("Instrument/EmpiricalOrthogonalFunction", 1 - 1),
               order = 2,
               creator = ConfigCommon.empirical_orthogonal_function,
               retrieve_bands = { true, true, true },
            },
            radiance_scaling = {
               apriori = ConfigCommon.hdf_apriori_i("Instrument/RadianceScaling/Coxmunk"),
               covariance = ConfigCommon.hdf_covariance_i("Instrument/RadianceScaling/Coxmunk"),
               creator = ConfigCommon.radiance_scaling_sv_fit_coxmunk_only,
               retrieve_bands = { true, true, true },
            },
         },
      },
      spectrum_effect = {

         creator = AcosConfig.spectrum_effect_list_h_and_m,
         speceff_h_gain = { "solar_model", "fluorescence" },
	 speceff_m_gain = { "solar_model",  },
         solar_model = {
            creator = ConfigCommon.solar_absorption_and_continuum,
            doppler_shift = {
               creator = ConfigCommon.solar_doppler_from_l1b,
               do_doppler_shift = true,
            },
            solar_absorption = {
               creator = ConfigCommon.solar_absorption_table,
            },
            solar_continuum = {
               creator = ConfigCommon.solar_continuum_table,
               convert_from_photon = true,
            },
         },
         fluorescence = {
            apriori = ConfigCommon.hdf_apriori("Fluorescence"),
            covariance = ConfigCommon.fluorescence_covariance("Fluorescence"),
            creator = ConfigCommon.fluorescence_effect_lambertian_only,
            reference_point = ConfigCommon.hdf_read_double_with_unit("Fluorescence/reference_point"),
            retrieved = true,
         },
      },
      spec_samp = {
         creator = ConfigCommon.nonuniform_spectrum_sampling,
         high_resolution_spectrum_spacing = DoubleWithUnit(0.01, "cm^-1"),
         nonunif_rt_grid_files = { 
            o2 = ConfigCommon.hdf_read_spec_dom("Spectrum_Sampling/nonuniform_grid_1"),
            weak_co2 = ConfigCommon.hdf_read_spec_dom("Spectrum_Sampling/nonuniform_grid_2"),
            strong_co2 = ConfigCommon.hdf_read_spec_dom("Spectrum_Sampling/nonuniform_grid_3"),
         },
      },
      rt = {
         creator = ConfigCommon.radiative_transfer_lsi,
         nadir_threshold = 0.25,
         lsi_constant = {
            dedicated_twostream = true,
            low_stream = 1, 
            high_stream = 8
         },
      },
      state_vector = {
         creator = ConfigCommon.state_vector_creator,
      },
      atmosphere = {
         creator = ConfigCommon.atmosphere_oco,
         constants = {
            creator = ConfigCommon.default_constant,
         },
         pressure = {
            apriori = ConfigCommon.ecmwf_pressure,
            covariance = ConfigCommon.hdf_covariance("Surface_Pressure"),
            a = ConfigCommon.hdf_read_double_1d("Pressure/Pressure_sigma_a"),
            b = ConfigCommon.hdf_read_double_1d("Pressure/Pressure_sigma_b"),
            creator = ConfigCommon.pressure_sigma,
         },
         temperature = {
            apriori = ConfigCommon.hdf_apriori("Temperature/Offset"),
            covariance = ConfigCommon.hdf_covariance("Temperature/Offset"),
            creator = ConfigCommon.temperature_ecmwf,
         },
         ground = {
            -- Pure lambertian
            lambertian = {
               apriori = AcosConfig.gosat_albedo_from_radiance(1),  
               covariance = ConfigCommon.hdf_covariance_i("Ground/Albedo"),
               retrieve_bands = { true, true, true },
               creator = ConfigCommon.lambertian_retrieval,
            },

            -- Coxmunk windspeed and refractive index inputs
            coxmunk = {
               refractive_index = ConfigCommon.hdf_apriori("Ground/Refractive_Index"),
               apriori = ConfigCommon.ecmwf_windspeed,
               covariance = ConfigCommon.hdf_covariance("Ground/Windspeed"),
               creator = ConfigCommon.coxmunk_retrieval,
            },

            -- Lambertian component of coxmunk + lambertian
            coxmunk_lambertian = {
               apriori = ConfigCommon.hdf_apriori_i("Ground/Coxmunk_Albedo"),
               covariance = ConfigCommon.hdf_covariance_i("Ground/Coxmunk_Albedo"),
               retrieve_bands = { true, true, true },
               creator = ConfigCommon.lambertian_retrieval,
            },
            
            -- Breon vegetative kernel with Rahman retrieved parameters
            breon_veg = {
               apriori = ConfigCommon.hdf_apriori("Ground/Rahman"),
               covariance = ConfigCommon.hdf_covariance("Ground/Rahman"),
               creator = ConfigCommon.breon_veg_retrieval,
            },
            
            -- Breon soil kernel with Rahman retrieved parameters
            breon_soil = {
               apriori = ConfigCommon.hdf_apriori("Ground/Rahman"),
               covariance = ConfigCommon.hdf_covariance("Ground/Rahman"),
               creator = ConfigCommon.breon_soil_retrieval,
            },

            creator = AcosConfig.ground_land_fraction,
         },
         aerosol = {
	    creator = ConfigCommon.merra_aerosol_creator,
	    max_aod = 0.2,
	    exp_aod = 0.8,
	    min_types = 2,
	    max_types = 2,
	    max_residual = 0.005,
	    apriori = ConfigCommon.hdf_apriori("/Aerosol/Merra/Gaussian/Log"),
	    covariance = ConfigCommon.hdf_covariance("/Aerosol/Merra/Gaussian/Log"),
	    -- Lua doesn't preserve order in a table, so we have a list
	    -- saying what order we want the Aerosols in
	    aerosols = {"Ice", "Water"},
	    Water = {
	       creator = ConfigCommon.aerosol_log_shape_gaussian,
	       apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
	       covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
	       property = ConfigCommon.hdf_aerosol_property("wc_008"),
	    },
	    Ice = {
	       creator = ConfigCommon.aerosol_log_shape_gaussian,
	       apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol", "Gaussian/Log"),
	       covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol", "Gaussian/Log"),
	       property = ConfigCommon.hdf_aerosol_property("ice_cloud_MODIS6_deltaM_1000"),
	    },
         },
         absorber = {
            creator = ConfigCommon.absorber_creator,
            gases = {"CO2", "H2O", "O2"},
            CO2 = {
               apriori = ConfigCommon.tccon_co2_apriori_ecmwf,
               covariance = ConfigCommon.hdf_covariance("Gas/CO2"),
               absco = "v4.2.0_unscaled/co2_v4.2.0_with_ctm.hdf",
	       table_scale = {1.0, 1.0038, 0.9946},
               creator = ConfigCommon.vmr_level,
            },
            H2O = {
               scale_apriori = 1.0,
               scale_cov = 0.25,
               absco = "v4.2.0_unscaled/h2o_v4.2.0.hdf",
               creator = ConfigCommon.vmr_ecmwf,
            },
            O2 = {
               apriori = ConfigCommon.hdf_read_double_1d("Gas/O2/average_mole_fraction"),
               absco = "v4.2.0_unscaled/o2_v4.2.0_drouin.hdf",
	       table_scale = 1.0125,
               creator = ConfigCommon.vmr_level_constant_well_mixed,
            },
         },
         altitude = {
            creator = ConfigCommon.hydrostatic_altitude,
         },
      },
   },
}