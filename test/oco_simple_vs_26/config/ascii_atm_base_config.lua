-- This sets up the "standard" run we use in the unit tests, you can then
-- use this and override any feature you want to be done differently
--
-- Pull anything from CommonConfig through the CommonConfig name and
-- anything in OcoConfig from that name instead of all through OcoConfig
-- to better differentiate where routines are being loaded from.

require "oco_config"

function albedo_ascii_file(self, i)
     return ConfigCommon.ascii_ground_apriori(self.config.albedo_file)(self, i)
end

function old_oco_psurf_cov_value(self) 
    local cov = Blitz_double_array_2d(1,1)
    cov:set(0, 0, 625)
    return cov
end

AsciiAtmBaseConfig = OcoConfig:new {
------------------------------------------------------------
--- Various constants used to describe input data.
------------------------------------------------------------

   -- As a convenience, we pass in a few fields by environment
   -- variables.
   sid_string = os.getenv("sounding_id"),
   spectrum_file = os.getenv("spectrum_file"),
   atmosphere_file = os.getenv("atmosphere_file"),
   albedo_file = os.getenv("albedo_file"),

   static_file = "../input/l2_oco_static_input.h5",
   static_solar_file = "../input/l2_solar_model.h5",

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

   solver = { max_iteration=12, max_divergence=8,
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
      l1b = {
         creator = OcoConfig.level1b_hdf,
         noise = {
            creator = OcoConfig.oco_noise,
            max_ms = { 7.00e20, 2.45e20, 1.25e20 },
         },
      },
      stokes_coefficient = {
	 creator = ConfigCommon.stokes_coefficient_constant,
	 value = ConfigCommon.stokes_coefficient_l1b,
      },
      instrument = {
         creator = ConfigCommon.ils_instrument,
         ils_half_width = { DoubleWithUnit(1.2e-3, "um"), 
                            DoubleWithUnit(5.1e-3, "um"),
                            DoubleWithUnit(8.4e-3, "um") },
         dispersion = {
            creator = ConfigCommon.dispersion_polynomial,
            apriori = ConfigCommon.l1b_spectral_coefficient_i,
            covariance = ConfigCommon.hdf_covariance_i("Instrument/Dispersion"),
            number_pixel = ConfigCommon.hdf_read_int_1d("Instrument/Dispersion/number_pixel"),
            retrieved = false,
            is_one_based = true,
         },
         ils_func = {
            creator = OcoConfig.ils_table_l1b,
         },
         instrument_correction = {
            creator = ConfigCommon.instrument_correction_list,
            ic = { },
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
         creator = ConfigCommon.uniform_spectrum_sampling,
         high_resolution_spectrum_spacing = DoubleWithUnit(0.01, "cm^-1"),
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
            pressure_levels = ConfigCommon.read_atmosphere_file("Pressure"),
            apriori = ConfigCommon.surface_pressure_from_atmosphere_file(),
            covariance = old_oco_psurf_cov_value, 
            creator = ConfigCommon.pressure_fixed_level,
         },
         temperature = {
            creator = ConfigCommon.temperature_fixed_level,
            levels = { apriori=ConfigCommon.read_atmosphere_file("T"), },
            offset = { apriori =ConfigCommon.hdf_apriori("Temperature/Offset"),
                       covariance = ConfigCommon.hdf_covariance("Temperature/Offset"), },
         },
         ground = {
            -- Pure lambertian
            lambertian = {
               apriori = albedo_ascii_file,
               covariance = ConfigCommon.hdf_covariance_i("Ground/Albedo"),
               retrieve_bands = { true, true, true },
               creator = ConfigCommon.lambertian_retrieval,
            },

            -- Coxmunk windspeed and refractive index inputs
            coxmunk = {
               refractive_index = ConfigCommon.hdf_apriori("Ground/Refractive_Index"),
               apriori = ConfigCommon.met_windspeed, 
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

            creator = OcoConfig.ground_land_water_indicator,
         },
         aerosol = {
            -- Add reading ascii aerosol later
            creator = ConfigCommon.rayleigh_only,
         },
         absorber = {
            creator = ConfigCommon.absorber_creator,
            gases = {"CO2", "H2O", "O2"},
            CO2 = {
               apriori = ConfigCommon.read_atmosphere_file("CO2"),
               covariance = ConfigCommon.hdf_covariance("Gas/CO2"),
               absco = "v3.3.0/lowres/co2_v3.3.0-lowres.hdf",
               creator = ConfigCommon.vmr_level,
            },
            H2O = {
               apriori = ConfigCommon.read_atmosphere_file("H2O"),
               scale_apriori = 1.0,
               scale_cov = 0.25,
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
