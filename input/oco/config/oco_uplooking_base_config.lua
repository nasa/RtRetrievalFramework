-- This sets up the "standard" run we use in the unit tests, you can then
-- use this and override any feature you want to be done differently
--
-- Pull anything from CommonConfig through the CommonConfig name and
-- anything in OcoConfig from that name instead of all through OcoConfig
-- to better differentiate where routines are being loaded from.

require "oco_config"
oco_base_config_dir = ConfigCommon.local_dir()

OcoUplookingBaseConfig = OcoConfig:new {

------------------------------------------------------------
--- Various constants used to describe input data.
------------------------------------------------------------

   -- As a convenience, we pass in a few fields by environment
   -- variables.
   sid_string = os.getenv("sounding_id"),
   spectrum_file = os.getenv("spectrum_file"),
   atmosphere_file = os.getenv("atmosphere_file"),

   static_file = oco_base_config_dir .. "/../input/l2_oco_static_input.h5",
   static_solar_file = config_common_dir .. "/../input/l2_solar_model.h5",

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

   solver = { max_iteration=10, max_divergence=2,
              max_chisq=1.4, threshold=0.1, gamma_initial=5.0,
              create = ConfigCommon.connor_solver, },

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
      input = {
          creator = ConfigCommon.l1b_input,
          l1b = {
             creator = OcoConfig.level1b_hdf,
             noise = {
                creator = OcoConfig.oco_noise,
                max_ms = { 7.00e20, 2.45e20, 1.25e20 },
             },
          },
      },
      stokes_coefficient = {
         creator = ConfigCommon.stokes_coefficient_constant,
         value = ConfigCommon.stokes_coefficient_l1b,
      },
      instrument = {
         creator = ConfigCommon.ils_instrument,
         ils_half_width = { DoubleWithUnit(4.09e-04, "um"), 
                            DoubleWithUnit(1.08e-03, "um"),
                            DoubleWithUnit(1.40e-03, "um") },
         dispersion = {
            creator = ConfigCommon.dispersion_polynomial,
            apriori = ConfigCommon.l1b_spectral_coefficient_i,
            covariance = ConfigCommon.hdf_covariance_i("Instrument/Dispersion"),
            number_pixel = ConfigCommon.hdf_read_int_1d("Instrument/Dispersion/number_pixel"),
            retrieved = true,
        is_one_based = true,
         },
         ils_func = {
            creator = OcoConfig.ils_table_l1b,
         },
         instrument_correction = {
            creator = ConfigCommon.instrument_correction_list,
            ic = { "radiance_scaling" },
            radiance_scaling = {
               apriori = ConfigCommon.hdf_apriori_i("Instrument/RadianceScaling/Uplooking"),
               covariance = ConfigCommon.hdf_covariance_i("Instrument/RadianceScaling/Uplooking"),
               creator = ConfigCommon.radiance_scaling_sv_fit,
               retrieve_bands = { true, true, true },
            },
         },
      },
      spectrum_effect = {
                creator = ConfigCommon.spectrum_effect_list,
                speceff = { "solar_model", },
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
               convert_from_photon = false,
            },
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
         creator = ConfigCommon.chapman_boa_rt,
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
            pressure_levels = ConfigCommon.read_atmosphere_file("Pressure"),
            apriori = ConfigCommon.surface_pressure_from_atmosphere_file(),
            covariance = ConfigCommon.hdf_covariance("Surface_Pressure"),
            creator = ConfigCommon.pressure_sigma_profile,
         },
         temperature = {
            temperature_levels = ConfigCommon.read_atmosphere_file("T"),
            apriori = ConfigCommon.hdf_apriori("Temperature/Offset"),
            covariance = ConfigCommon.hdf_covariance("Temperature/Offset"),
            creator = ConfigCommon.temperature_level_offset,
         },
         ground = {
            creator = ConfigCommon.no_ground,
         },
         aerosol = {
            -- No aerosols for uplooking mode
            creator = ConfigCommon.rayleigh_only,
         },
         absorber = {
            creator = ConfigCommon.absorber_creator,
            gases = {"CO2", "H2O", "O2"},
            CO2 = {
               apriori = ConfigCommon.co2_from_atmosphere_or_tccon("CO2"),
               covariance = ConfigCommon.hdf_covariance("Gas/CO2"),
               absco = "v5.2_final/co2_v52.hdf",
               table_scale = {1.0, 1.0, 1.004},
               creator = ConfigCommon.vmr_level,
            },
            H2O = {
               vmr_profile = ConfigCommon.read_atmosphere_file("H2O"),
               scale_apriori = 1.0,
               scale_cov = 0.1,
               absco = "v5.2_final/h2o_v52.hdf",
               creator = ConfigCommon.vmr_level_scaled,
            },
            O2 = {
               apriori = ConfigCommon.hdf_read_double_1d("Gas/O2/average_mole_fraction"),
               absco = "v5.2_final/o2_v52.hdf",
               table_scale = 1.0048,
               creator = ConfigCommon.vmr_level_constant_well_mixed,
            },
         },
         altitude = {
            creator = ConfigCommon.hydrostatic_altitude,
         },
	 relative_humidity = {
	    creator = ConfigCommon.calc_relative_humidity,
	 },
      },
   },
}
