-- This sets up the "standard" run we use in the unit tests, you can then
-- use this and override any feature you want to be done differently
--
-- Pull anything from CommonConfig through the CommonConfig name and
-- anything in FtsConfig from that name instead of all through FtsConfig
-- to better differentiate where routines are being loaded from.

require "fts_unit_test_config"
fts_test_base_config_dir = ConfigCommon.local_dir()

FtsTestBaseConfig = FtsUnitTestConfig:new {
------------------------------------------------------------
--- Various constants used to describe input data.
------------------------------------------------------------

   -- As a convenience, we pass in a few fields by environment
   -- variables.
   runlog_file = fts_test_base_config_dir .. "/in/tccon_runlog.grl",
   static_file = fts_test_base_config_dir .. 
      "/l2_fts_unit_test_static_input.h5",
   atmosphere_file = fts_test_base_config_dir .. 
      "/in/scene/atmosphere/atmosphere_fts.dat",

   -- For FTS runs, each spectrometer can have a different
   -- source L1B file
   spectrum_1_file = fts_test_base_config_dir .. 
      "/in/l1b/spec/pa20091103saaaaa_100223160344.008",

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

   absco_path        = "/groups/algorithm/l2_fp/absco",

------------------------------------------------------------
-- Connor solver
------------------------------------------------------------

   solver = { max_iteration=10, max_divergence=10,
	      max_chisq=1.0, threshold=0.1, gamma_initial=0.0,
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
	 creator = FtsUnitTestConfig.fts_spectral_window,
	 spec_win_file = "in/oco_l2.win",
      },
      input = {
         creator = ConfigCommon.l1b_input,
         l1b = {
   	 creator = FtsUnitTestConfig.level1b_fts,
         },
      },
      stokes_coefficient = {
	 creator = ConfigCommon.stokes_coefficient_constant,
	 value = ConfigCommon.stokes_coefficient_l1b,
      },
      instrument = {
	 creator = FtsUnitTestConfig.fts_instrument,
	 dispersion = {
	    creator = ConfigCommon.dispersion_polynomial,
	    apriori = ConfigCommon.l1b_spectral_coefficient_i,
	    covariance = ConfigCommon.hdf_covariance_i("Instrument/Dispersion"),
	    perturb = ConfigCommon.hdf_read_double_2d("Instrument/Dispersion/perturb"),
	    number_pixel = FtsUnitTestConfig.l1b_number_pixel,
	    is_one_based = false,
	 },
	 instrument_correction = {
	    creator = ConfigCommon.instrument_correction_list,
	    ic = { "continuum", },
	    continuum = {
	       creator = ConfigCommon.radiance_scaling_sv_fit,
	       apriori = ConfigCommon.hdf_apriori_i("Instrument/Continuum"),
	       covariance = ConfigCommon.hdf_covariance_i("Instrument/Continuum"),
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
	       creator = ConfigCommon.solar_doppler_from_runlog,
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
      },
      spec_samp = {
	 creator = ConfigCommon.uniform_spectrum_sampling,
	 high_resolution_spectrum_spacing = { DoubleWithUnit(0.01, "cm^-1"),
					      DoubleWithUnit(0.01, "cm^-1"),
					      DoubleWithUnit(0.01, "cm^-1") },
      },
      rt = {
	 creator = FtsUnitTestConfig.chapman_boa_rt_sza_calculate,
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
	    apriori = FtsUnitTestConfig.surface_pressure_from_runlog,
	    covariance = ConfigCommon.hdf_covariance("Surface_Pressure"),
	    creator = ConfigCommon.pressure_fixed_level,
	 },
	 temperature = {
	    levels = {
	       apriori = ConfigCommon.read_atmosphere_file("T"),
	       covariance = ConfigCommon.hdf_apriori("Temperature/Levels"),
	    },
	    offset = {
	       apriori = ConfigCommon.hdf_apriori("Temperature/Offset"),
	       covariance = ConfigCommon.hdf_covariance("Temperature/Offset"),
	    },
	    creator = ConfigCommon.temperature_fixed_level,
	 },
	 ground = {
	    creator = ConfigCommon.no_ground,
	 },
	 aerosol = {
	    -- No aerosols for FTS mode
	    creator = ConfigCommon.rayleigh_only,
	 },
	 absorber = {
	    creator = ConfigCommon.absorber_creator,
	    gases = {"CO2", "H2O", "O2"},
	    CO2 = {
	       apriori = ConfigCommon.read_atmosphere_file("CO2"),
	       scale_apriori = 1.0,
	       scale_cov = 0.1,
	       absco = "v3.3.0/lowres/co2_v3.3.0-lowres.hdf",
	       creator = ConfigCommon.vmr_fixed_level_scaled,
	    },
	    H2O = {
	       apriori = ConfigCommon.read_atmosphere_file("H2O"),
	       scale_apriori = 1.0,
	       scale_cov = 0.1,
	       absco = "v3.3.0/lowres/h2o_v3.3.0-lowres.hdf",
	       creator = ConfigCommon.vmr_fixed_level_scaled,
	    },
	    O2 = {
	       apriori = ConfigCommon.hdf_read_double_1d("Gas/O2/average_mole_fraction"),
	       absco = "v3.3.0/lowres/o2_v3.3.0-lowres.hdf",
	       creator = ConfigCommon.vmr_fixed_level_constant_well_mixed,
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
