-- This sets up the "standard" run we use in the unit tests, you can then
-- use this and override any feature you want to be done differently

--- This is found in the input/gosat/config directory
require "acos_config"
fixed_level_base_config_dir = ConfigCommon.local_dir()

FixedLevelBaseConfig = AcosConfig:new {
------------------------------------------------------------
--- Various constants used to describe input data.
------------------------------------------------------------

   sid_string = "20090627205850",
   cloud_file = fixed_level_base_config_dir .. "/in/cloud.h5",
   static_file = fixed_level_base_config_dir .. 
      "/l2_fixed_level_static_input.h5",
   soundinginfo_file = fixed_level_base_config_dir .. 
      "/in/l1b/soundinginfo.dat",
   spectrum_file     = fixed_level_base_config_dir ..
      "/in/l1b/spec/spectra.dat",
   absco_path        = "/groups/algorithm/l2_fp/absco",

------------------------------------------------------------
-- Set this true to get diagnostic messages to help debug
-- problems with Lua
------------------------------------------------------------

   diagnostic = false,

------------------------------------------------------------
-- Connor solver
------------------------------------------------------------

   solver = { threshold=1.0,
              min_iteration=1,
              max_iteration=12,
              max_divergence=8,
              max_chisq=1.4,
              gamma_initial=10.0,
              h2o_scale_index0=-20,
              h2o_scale_index1=-20,
              h2o_scale_cov_initial=0.001,
              ch4_scale_index0=-21,
              ch4_scale_index1=-21,
              ch4_scale_cov_initial=0.001,
              co_scale_index0=-22,
              co_scale_index1=-22,
              co_scale_cov_initial=0.0001,
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
         creator = ConfigCommon.spectral_window_hdf,
      },
      instrument = {
         creator = ConfigCommon.ils_instrument,
         ils_half_width = { DoubleWithUnit(20, "cm^-1"),
                            DoubleWithUnit(20, "cm^-1"),
                            DoubleWithUnit(20, "cm^-1") },
         ils_func = {
            creator = ConfigCommon.ils_table,
         },
         dispersion = {
            creator = ConfigCommon.dispersion_polynomial,
            apriori = ConfigCommon.hdf_apriori_with_unit_i("Instrument/Dispersion"),
            covariance = ConfigCommon.hdf_covariance_i("Instrument/Dispersion"),
            number_pixel = ConfigCommon.hdf_read_int_1d("Instrument/Dispersion/number_pixel"),
            is_one_based = true,
         },
         instrument_correction = {
            creator = ConfigCommon.instrument_correction_list,
            ic = {},
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
               convert_from_photon = true,
            },
         },
      },
      input = {
         creator = ConfigCommon.l1b_input,
         l1b = {
            creator = ConfigCommon.level1b_ascii,
            noise = {
               creator = ConfigCommon.noise_ascii,
            },
         },
      },
      stokes_coefficient = {
	 creator = ConfigCommon.stokes_coefficient_constant,
	 value = ConfigCommon.stokes_coefficient_l1b,
      },
      spec_samp = {
         creator = ConfigCommon.uniform_spectrum_sampling,
         high_resolution_spectrum_spacing = DoubleWithUnit(0.01, "cm^-1"),
         nonunif_rt_grid_files = { 
            o2 = "common/nonunif_rt_grid__gosat_abo2_oco__absco_v3.1.0__wn.dat",
            weak_co2 = "common/nonunif_rt_grid__gosat_wco2_oco__absco_v3.1.0__wn.dat",
            strong_co2 = "common/nonunif_rt_grid__gosat_sco2_oco__absco_v3.1.0__wn.dat",
         },
      },
      rt = {
         creator = ConfigCommon.radiative_transfer_lsi,
         nadir_threshold = 0.25,
         lsi_constant = {
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
            pressure_levels = ConfigCommon.hdf_read_double_1d("Pressure/Pressure"),
            apriori = ConfigCommon.hdf_apriori("Surface_Pressure"),
            covariance = ConfigCommon.hdf_covariance("Surface_Pressure"),
            creator = ConfigCommon.pressure_fixed_level,
         },
         temperature = {
            levels = {
               apriori = ConfigCommon.hdf_apriori("Temperature/Levels"),
               covariance = ConfigCommon.hdf_apriori("Temperature/Levels"),
            },
            offset = {
               apriori = ConfigCommon.hdf_apriori("Temperature/Offset"),
               covariance = ConfigCommon.hdf_covariance("Temperature/Offset"),
            },
            creator = ConfigCommon.temperature_fixed_level,
         },
         ground = {
            -- Pure lambertian
            lambertian = {
               apriori = ConfigCommon.hdf_apriori_i("Ground/Albedo"),
               covariance = ConfigCommon.hdf_covariance_i("Ground/Albedo"),
               retrieve_bands = { true, true, true },
               creator = ConfigCommon.lambertian_retrieval,
            },

            -- Coxmunk windspeed and refractive index inputs
            coxmunk = {
               refractive_index = ConfigCommon.hdf_apriori("Ground/Refractive_Index"),
               apriori = ConfigCommon.hdf_apriori("Ground/Windspeed"),
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

            -- Lambertian component of coxmunk + lambertian
            coxmunk_scaled = {
               apriori = ConfigCommon.hdf_apriori_i("Ground/Coxmunk_Scaled"),
               covariance = ConfigCommon.hdf_covariance_i("Ground/Coxmunk_Scaled"),
               retrieve_bands = { true, true, true },
               scaled_brdf_name = "CoxMunk",
               creator = ConfigCommon.brdf_scale_retrieval,
            },

            -- Brdf vegetative kernel with Rahman retrieved parameters
            brdf_veg = {
               apriori = ConfigCommon.hdf_apriori_i("Ground/Brdf"),
               covariance = ConfigCommon.hdf_covariance_i("Ground/Brdf"),
               creator = ConfigCommon.brdf_veg_retrieval,
            },
            
            -- Brdf soil kernel with Rahman retrieved parameters
            brdf_soil = {
               apriori = ConfigCommon.hdf_apriori_i("Ground/Brdf"),
               covariance = ConfigCommon.hdf_covariance_i("Ground/Brdf"),
               creator = ConfigCommon.brdf_soil_retrieval,
            },

            creator = ConfigCommon.ground_lambertian,
         },
         aerosol = {
            creator = ConfigCommon.aerosol_creator,
            -- Lua doesn't preserve order in a table, so we have a list
            -- saying what order we want the Aerosols in
            aerosols = {"Kahn_2b", "Kahn_3b", "Water", "Ice"},
            Kahn_2b = {
               creator = ConfigCommon.aerosol_log_profile,
               apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol"),
               covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol"),
               property = ConfigCommon.hdf_aerosol_property("Aerosol/Kahn_2b"),
            },
            Kahn_3b = {
               creator = ConfigCommon.aerosol_log_profile,
               apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol"),
               covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol"),
               property = ConfigCommon.hdf_aerosol_property("Aerosol/Kahn_3b"),
            },
            Water = {
               creator = ConfigCommon.aerosol_log_profile,
               apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol"),
               covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol"),
               property = ConfigCommon.hdf_aerosol_property("Aerosol/Water"),
            },
            Ice = {
               creator = ConfigCommon.aerosol_log_profile,
               apriori = ConfigCommon.hdf_aerosol_apriori("Aerosol"),
               covariance = ConfigCommon.hdf_aerosol_covariance("Aerosol"),
               property = ConfigCommon.hdf_aerosol_property("Aerosol/Ice"),
            },
         },
         absorber = {
            creator = ConfigCommon.absorber_creator,
            gases = {"CO2", "H2O", "O2"},
            CO2 = {
               apriori = ConfigCommon.hdf_apriori("Gas/CO2"),
               covariance = ConfigCommon.hdf_covariance("Gas/CO2"),
               absco = "v3.3.0/lowres/co2_v3.3.0-lowres.hdf",
               creator = ConfigCommon.vmr_fixed_level,
            },
            H2O = {
               apriori = ConfigCommon.hdf_apriori("Gas/H2O"),
               scale_apriori = 1.0,
               scale_cov = 0.25,
               absco = "v3.3.0/lowres/h2o_v3.3.0-lowres.hdf",
               creator = ConfigCommon.vmr_fixed_level_scaled,
            },
            O2 = {
               apriori = ConfigCommon.hdf_apriori("Gas/O2"),
               absco = "v3.3.0/lowres/o2_v3.3.0-lowres.hdf",
               creator = ConfigCommon.vmr_fixed_level_constant,
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
