#include <blitz/array.h>
#include <cmath>
#include "unit_test_support.h"
#include "lidort_interface_types.h"
#include "spurr_brdf_types.h"

// Default maximum number of layers
#define LIDORT_MAXLAYER 23

using namespace FullPhysics;
using namespace blitz;

/* This file was auto-generated */

BOOST_FIXTURE_TEST_SUITE(lidort_interface_types, GlobalFixture)


BOOST_AUTO_TEST_CASE(lidort_pars)
{
  // Test obtaining instance
  Lidort_Pars tst_obj = Lidort_Pars::instance();
  // Test that the value obtained from fortran is as expected
  BOOST_CHECK_EQUAL(tst_obj.lidort_version_number, "3.6");
  BOOST_CHECK_EQUAL(tst_obj.lidort_inunit, 21);
  BOOST_CHECK_EQUAL(tst_obj.lidort_scenunit, 22);
  BOOST_CHECK_EQUAL(tst_obj.lidort_funit, 23);
  BOOST_CHECK_EQUAL(tst_obj.lidort_resunit, 24);
  BOOST_CHECK_EQUAL(tst_obj.lidort_errunit, 25);
  BOOST_CHECK_EQUAL(tst_obj.lidort_dbgunit, 71);
  BOOST_CHECK_EQUAL(tst_obj.max_messages, 25);
  BOOST_CHECK_EQUAL(tst_obj.maxthreads, 8);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams, 10);
  BOOST_CHECK_EQUAL(tst_obj.maxlayers, LIDORT_MAXLAYER);
  BOOST_CHECK_EQUAL(tst_obj.maxfinelayers, 4);
  BOOST_CHECK_EQUAL(tst_obj.maxmoments_input, 180);
  BOOST_CHECK_EQUAL(tst_obj.max_thermal_coeffs, 3);
  BOOST_CHECK_EQUAL(tst_obj.maxbeams, 4);
  BOOST_CHECK_EQUAL(tst_obj.max_user_streams, 4);
  BOOST_CHECK_EQUAL(tst_obj.max_user_relazms, 3);
  BOOST_CHECK_EQUAL(tst_obj.max_user_levels, 5);
  BOOST_CHECK_EQUAL(tst_obj.max_partlayers, 2);
  BOOST_CHECK_EQUAL(tst_obj.max_directions, 2);
  BOOST_CHECK_EQUAL(tst_obj.max_brdf_kernels, 3);
  BOOST_CHECK_EQUAL(tst_obj.max_brdf_parameters, 3);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams_brdf, 101);
  BOOST_CHECK_EQUAL(tst_obj.max_msrs_muquad, 50);
  BOOST_CHECK_EQUAL(tst_obj.max_msrs_phiquad, 100);
  BOOST_CHECK_EQUAL(tst_obj.max_atmoswfs, 7);
  BOOST_CHECK_EQUAL(tst_obj.max_surfacewfs, 7);
  BOOST_CHECK_EQUAL(tst_obj.max_sleavewfs, 2);
  BOOST_CHECK_EQUAL(tst_obj.max_geometries, tst_obj.max_user_streams*tst_obj.max_user_relazms*tst_obj.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.max_allstrms, tst_obj.max_user_streams + tst_obj.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.max_allstrms_p1, tst_obj.max_allstrms + tst_obj.maxbeams*tst_obj.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.maxmoments, 2*tst_obj.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.maxfourier, 2*tst_obj.maxstreams - 1);
  BOOST_CHECK_EQUAL(tst_obj.maxsthalf_brdf, tst_obj.maxstreams_brdf / 2);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams_2, 2*tst_obj.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.maxstreams_p1, tst_obj.maxstreams + 1);
  BOOST_CHECK_EQUAL(tst_obj.maxtotal, tst_obj.maxlayers*tst_obj.maxstreams_2);
  BOOST_CHECK_EQUAL(tst_obj.maxbandtotal, 9*tst_obj.maxstreams - 2);
  BOOST_CHECK_CLOSE(tst_obj.one, 1.0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.zero, 0.0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.onep5, 1.5, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.two, 2.0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.three, 3.0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.four, 4.0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.quarter, 0.25, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.half, 0.5, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.minus_one, -tst_obj.one, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.minus_two, -tst_obj.two, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.pie, 3.141592653589793e0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.deg_to_rad, tst_obj.pie/180.0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.pi2, tst_obj.two  * tst_obj.pie, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.pi4, tst_obj.four * tst_obj.pie, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.pio2, tst_obj.half * tst_obj.pie, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.pio4, tst_obj.quarter * tst_obj.pie, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.eps3, 0.001, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.eps4, 0.0001, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.eps5, 0.00001, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.smallnum, 0.000000001, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.bigexp, 32.0, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.hopital_tolerance, tst_obj.eps5, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.omega_smallnum, 0.00000001, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.max_tau_spath, tst_obj.bigexp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.max_tau_upath, tst_obj.bigexp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.max_tau_qpath, tst_obj.bigexp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.lidort_serious, 4);
  BOOST_CHECK_EQUAL(tst_obj.lidort_warning, 3);
  BOOST_CHECK_EQUAL(tst_obj.lidort_info, 2);
  BOOST_CHECK_EQUAL(tst_obj.lidort_debug, 1);
  BOOST_CHECK_EQUAL(tst_obj.lidort_success, 0);
  BOOST_CHECK_EQUAL(tst_obj.upidx, 1);
  BOOST_CHECK_EQUAL(tst_obj.dnidx, 2);
  BOOST_CHECK_EQUAL(tst_obj.lambertian_idx, LAMBERTIAN);
  BOOST_CHECK_EQUAL(tst_obj.rossthin_idx, ROSSTHIN);
  BOOST_CHECK_EQUAL(tst_obj.rossthick_idx, ROSSTHICK);
  BOOST_CHECK_EQUAL(tst_obj.lisparse_idx, LISPARSE);
  BOOST_CHECK_EQUAL(tst_obj.lidense_idx, LIDENSE);
  BOOST_CHECK_EQUAL(tst_obj.hapke_idx, HAPKE);
  BOOST_CHECK_EQUAL(tst_obj.roujean_idx, ROUJEAN);
  BOOST_CHECK_EQUAL(tst_obj.rahman_idx, RAHMAN);
  BOOST_CHECK_EQUAL(tst_obj.coxmunk_idx, COXMUNK);
  BOOST_CHECK_EQUAL(tst_obj.breonveg_idx, BREONVEG);
  BOOST_CHECK_EQUAL(tst_obj.breonsoil_idx, BREONSOIL);
  BOOST_CHECK_EQUAL(tst_obj.maxbrdf_idx, tst_obj.breonsoil_idx);
  
}

BOOST_AUTO_TEST_CASE(brdf_linsup_inputs)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Brdf_Linsup_Inputs tst_obj = Brdf_Linsup_Inputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.bs_do_kernel_factor_wfs().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_kernel_params_wfs().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_kernel_params_wfs().extent(1), lid_pars.max_brdf_parameters);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_kparams_derivs().extent(0), lid_pars.max_brdf_kernels);
  

  // Test initialization
  blitz::Array<bool, 1> bs_do_kernel_factor_wfs_exp(tst_obj.bs_do_kernel_factor_wfs().shape());
  bs_do_kernel_factor_wfs_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_do_kernel_factor_wfs(), bs_do_kernel_factor_wfs_exp, 1e-10);
  blitz::Array<bool, 2> bs_do_kernel_params_wfs_exp(tst_obj.bs_do_kernel_params_wfs().shape());
  bs_do_kernel_params_wfs_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_do_kernel_params_wfs(), bs_do_kernel_params_wfs_exp, 1e-10);
  blitz::Array<bool, 1> bs_do_kparams_derivs_exp(tst_obj.bs_do_kparams_derivs().shape());
  bs_do_kparams_derivs_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_do_kparams_derivs(), bs_do_kparams_derivs_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_surface_wfs(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_kernel_factor_wfs(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_kernel_params_wfs(), 0);
  
}

BOOST_AUTO_TEST_CASE(brdf_linsup_outputs)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Brdf_Linsup_Outputs tst_obj = Brdf_Linsup_Outputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_exactdb_brdfunc().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_exactdb_brdfunc().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_exactdb_brdfunc().extent(2), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_exactdb_brdfunc().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f_0().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f_0().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_brdf_f().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f_0().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f_0().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_brdf_f().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_emissivity().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_emissivity().extent(1), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_emissivity().extent(2), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_emissivity().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_emissivity().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_ls_user_emissivity().extent(2), lid_pars.maxthreads);
  

  // Test initialization
  blitz::Array<double, 4> bs_ls_exactdb_brdfunc_exp(tst_obj.bs_ls_exactdb_brdfunc().shape());
  bs_ls_exactdb_brdfunc_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_exactdb_brdfunc(), bs_ls_exactdb_brdfunc_exp, 1e-10);
  blitz::Array<double, 4> bs_ls_brdf_f_0_exp(tst_obj.bs_ls_brdf_f_0().shape());
  bs_ls_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_brdf_f_0(), bs_ls_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 4> bs_ls_brdf_f_exp(tst_obj.bs_ls_brdf_f().shape());
  bs_ls_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_brdf_f(), bs_ls_brdf_f_exp, 1e-10);
  blitz::Array<double, 4> bs_ls_user_brdf_f_0_exp(tst_obj.bs_ls_user_brdf_f_0().shape());
  bs_ls_user_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_user_brdf_f_0(), bs_ls_user_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 4> bs_ls_user_brdf_f_exp(tst_obj.bs_ls_user_brdf_f().shape());
  bs_ls_user_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_user_brdf_f(), bs_ls_user_brdf_f_exp, 1e-10);
  blitz::Array<double, 3> bs_ls_emissivity_exp(tst_obj.bs_ls_emissivity().shape());
  bs_ls_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_emissivity(), bs_ls_emissivity_exp, 1e-10);
  blitz::Array<double, 3> bs_ls_user_emissivity_exp(tst_obj.bs_ls_user_emissivity().shape());
  bs_ls_user_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_ls_user_emissivity(), bs_ls_user_emissivity_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(brdf_sup_inputs)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Brdf_Sup_Inputs tst_obj = Brdf_Sup_Inputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.bs_beam_szas().extent(0), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_relazms().extent(0), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_angles_input().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_which_brdf().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_brdf_parameters().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_parameters().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_parameters().extent(1), lid_pars.max_brdf_parameters);
  BOOST_CHECK_EQUAL(tst_obj.bs_lambertian_kernel_flag().extent(0), lid_pars.max_brdf_kernels);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_factors().extent(0), lid_pars.max_brdf_kernels);
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.bs_do_user_streams(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_brdf_surface(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_surface_emission(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_nstreams(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_nbeams(), 0);
  blitz::Array<double, 1> bs_beam_szas_exp(tst_obj.bs_beam_szas().shape());
  bs_beam_szas_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_beam_szas(), bs_beam_szas_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_user_relazms(), 0);
  blitz::Array<double, 1> bs_user_relazms_exp(tst_obj.bs_user_relazms().shape());
  bs_user_relazms_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_relazms(), bs_user_relazms_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_user_streams(), 0);
  blitz::Array<double, 1> bs_user_angles_input_exp(tst_obj.bs_user_angles_input().shape());
  bs_user_angles_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_angles_input(), bs_user_angles_input_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_n_brdf_kernels(), 0);
  blitz::Array<int, 1> bs_which_brdf_exp(tst_obj.bs_which_brdf().shape());
  bs_which_brdf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_which_brdf(), bs_which_brdf_exp, 1e-10);
  blitz::Array<int, 1> bs_n_brdf_parameters_exp(tst_obj.bs_n_brdf_parameters().shape());
  bs_n_brdf_parameters_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_n_brdf_parameters(), bs_n_brdf_parameters_exp, 1e-10);
  blitz::Array<double, 2> bs_brdf_parameters_exp(tst_obj.bs_brdf_parameters().shape());
  bs_brdf_parameters_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_brdf_parameters(), bs_brdf_parameters_exp, 1e-10);
  blitz::Array<bool, 1> bs_lambertian_kernel_flag_exp(tst_obj.bs_lambertian_kernel_flag().shape());
  bs_lambertian_kernel_flag_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_lambertian_kernel_flag(), bs_lambertian_kernel_flag_exp, 1e-10);
  blitz::Array<double, 1> bs_brdf_factors_exp(tst_obj.bs_brdf_factors().shape());
  bs_brdf_factors_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_brdf_factors(), bs_brdf_factors_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.bs_nstreams_brdf(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_shadow_effect(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_exactonly(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_glitter_msrcorr(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_do_glitter_msrcorr_exactonly(), false);
  BOOST_CHECK_EQUAL(tst_obj.bs_glitter_msrcorr_order(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_glitter_msrcorr_nmuquad(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_glitter_msrcorr_nphiquad(), 0);
  
}

BOOST_AUTO_TEST_CASE(brdf_sup_outputs)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Brdf_Sup_Outputs tst_obj = Brdf_Sup_Outputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.bs_exactdb_brdfunc().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_exactdb_brdfunc().extent(1), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.bs_exactdb_brdfunc().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f_0().extent(1), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f_0().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f().extent(1), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_brdf_f().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f_0().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f_0().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_brdf_f().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_emissivity().extent(0), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.bs_emissivity().extent(1), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_emissivity().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.bs_user_emissivity().extent(1), lid_pars.maxthreads);
  

  // Test initialization
  blitz::Array<double, 3> bs_exactdb_brdfunc_exp(tst_obj.bs_exactdb_brdfunc().shape());
  bs_exactdb_brdfunc_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_exactdb_brdfunc(), bs_exactdb_brdfunc_exp, 1e-10);
  blitz::Array<double, 3> bs_brdf_f_0_exp(tst_obj.bs_brdf_f_0().shape());
  bs_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_brdf_f_0(), bs_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 3> bs_brdf_f_exp(tst_obj.bs_brdf_f().shape());
  bs_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_brdf_f(), bs_brdf_f_exp, 1e-10);
  blitz::Array<double, 3> bs_user_brdf_f_0_exp(tst_obj.bs_user_brdf_f_0().shape());
  bs_user_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_brdf_f_0(), bs_user_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 3> bs_user_brdf_f_exp(tst_obj.bs_user_brdf_f().shape());
  bs_user_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_brdf_f(), bs_user_brdf_f_exp, 1e-10);
  blitz::Array<double, 2> bs_emissivity_exp(tst_obj.bs_emissivity().shape());
  bs_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_emissivity(), bs_emissivity_exp, 1e-10);
  blitz::Array<double, 2> bs_user_emissivity_exp(tst_obj.bs_user_emissivity().shape());
  bs_user_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.bs_user_emissivity(), bs_user_emissivity_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(brdf_input_exception_handling)
{
  // Test constructor
  Brdf_Input_Exception_Handling tst_obj = Brdf_Input_Exception_Handling();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.bs_status_inputread(), 0);
  BOOST_CHECK_EQUAL(tst_obj.bs_ninputmessages(), 0);
  
}

BOOST_AUTO_TEST_CASE(lidort_fixed_lincontrol)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Fixed_Lincontrol tst_obj = Lidort_Fixed_Lincontrol();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_layer_vary_flag().extent(0), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_layer_vary_number().extent(0), lid_pars.maxlayers);
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_do_column_linearization(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_profile_linearization(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_surface_linearization(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_sleave_wfs(), false);
  blitz::Array<bool, 1> ts_layer_vary_flag_exp(tst_obj.ts_layer_vary_flag().shape());
  ts_layer_vary_flag_exp = false;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_layer_vary_flag(), ts_layer_vary_flag_exp, 1e-10);
  blitz::Array<int, 1> ts_layer_vary_number_exp(tst_obj.ts_layer_vary_number().shape());
  ts_layer_vary_number_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_layer_vary_number(), ts_layer_vary_number_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_totalcolumn_wfs(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_surface_wfs(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_sleave_wfs(), 0);
  
}

BOOST_AUTO_TEST_CASE(lidort_fixed_linoptical)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Fixed_Linoptical tst_obj = Lidort_Fixed_Linoptical();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_l_deltau_vert_input().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_deltau_vert_input().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_deltau_vert_input().extent(2), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_omega_total_input().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_omega_total_input().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_omega_total_input().extent(2), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_phasmoms_total_input().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_phasmoms_total_input().extent(1), lid_pars.maxmoments_input+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_phasmoms_total_input().extent(2), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_l_phasmoms_total_input().extent(3), lid_pars.maxthreads);
  

  // Test initialization
  blitz::Array<double, 3> ts_l_deltau_vert_input_exp(tst_obj.ts_l_deltau_vert_input().shape());
  ts_l_deltau_vert_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_l_deltau_vert_input(), ts_l_deltau_vert_input_exp, 1e-10);
  blitz::Array<double, 3> ts_l_omega_total_input_exp(tst_obj.ts_l_omega_total_input().shape());
  ts_l_omega_total_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_l_omega_total_input(), ts_l_omega_total_input_exp, 1e-10);
  blitz::Array<double, 4> ts_l_phasmoms_total_input_exp(tst_obj.ts_l_phasmoms_total_input().shape());
  ts_l_phasmoms_total_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_l_phasmoms_total_input(), ts_l_phasmoms_total_input_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_fixed_lininputs)
{
  // Test constructor
  Lidort_Fixed_Lininputs tst_obj = Lidort_Fixed_Lininputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(lidort_modified_lininputs)
{
  // Test constructor
  Lidort_Modified_Lininputs tst_obj = Lidort_Modified_Lininputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.dummy(), 0);
  
}

BOOST_AUTO_TEST_CASE(lidort_linatmos)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Linatmos tst_obj = Lidort_Linatmos();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf().extent(4), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_columnwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_columnwf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_columnwf().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_columnwf().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_columnwf().extent(4), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_columnwf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_columnwf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_columnwf().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_columnwf().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_columnwf().extent(4), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(3), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf().extent(5), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_profilewf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_profilewf().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_profilewf().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_profilewf().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_profilewf().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_profilewf().extent(5), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_profilewf().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_profilewf().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_profilewf().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_profilewf().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_profilewf().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_profilewf().extent(5), lid_pars.maxthreads);
  

  // Test initialization
  blitz::Array<double, 5> ts_columnwf_exp(tst_obj.ts_columnwf().shape());
  ts_columnwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_columnwf(), ts_columnwf_exp, 1e-10);
  blitz::Array<double, 5> ts_mint_columnwf_exp(tst_obj.ts_mint_columnwf().shape());
  ts_mint_columnwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_mint_columnwf(), ts_mint_columnwf_exp, 1e-10);
  blitz::Array<double, 5> ts_flux_columnwf_exp(tst_obj.ts_flux_columnwf().shape());
  ts_flux_columnwf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_flux_columnwf(), ts_flux_columnwf_exp, 1e-10);
  blitz::Array<double, 6> ts_profilewf_exp(tst_obj.ts_profilewf().shape());
  ts_profilewf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_profilewf(), ts_profilewf_exp, 1e-10);
  blitz::Array<double, 6> ts_mint_profilewf_exp(tst_obj.ts_mint_profilewf().shape());
  ts_mint_profilewf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_mint_profilewf(), ts_mint_profilewf_exp, 1e-10);
  blitz::Array<double, 6> ts_flux_profilewf_exp(tst_obj.ts_flux_profilewf().shape());
  ts_flux_profilewf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_flux_profilewf(), ts_flux_profilewf_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_linsurf)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Linsurf tst_obj = Lidort_Linsurf();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf().extent(4), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_surfacewf().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_surfacewf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_surfacewf().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_surfacewf().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_mint_surfacewf().extent(4), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_surfacewf().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_surfacewf().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_surfacewf().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_surfacewf().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_surfacewf().extent(4), lid_pars.maxthreads);
  

  // Test initialization
  blitz::Array<double, 5> ts_surfacewf_exp(tst_obj.ts_surfacewf().shape());
  ts_surfacewf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_surfacewf(), ts_surfacewf_exp, 1e-10);
  blitz::Array<double, 5> ts_mint_surfacewf_exp(tst_obj.ts_mint_surfacewf().shape());
  ts_mint_surfacewf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_mint_surfacewf(), ts_mint_surfacewf_exp, 1e-10);
  blitz::Array<double, 5> ts_flux_surfacewf_exp(tst_obj.ts_flux_surfacewf().shape());
  ts_flux_surfacewf_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_flux_surfacewf(), ts_flux_surfacewf_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_linoutputs)
{
  // Test constructor
  Lidort_Linoutputs tst_obj = Lidort_Linoutputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(lidort_linsup_brdf)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Linsup_Brdf tst_obj = Lidort_Linsup_Brdf();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_exactdb_brdfunc().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_exactdb_brdfunc().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_exactdb_brdfunc().extent(2), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_exactdb_brdfunc().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f_0().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f_0().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_brdf_f().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f_0().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f_0().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_brdf_f().extent(3), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_emissivity().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_emissivity().extent(1), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_emissivity().extent(2), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_emissivity().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_emissivity().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_ls_user_emissivity().extent(2), lid_pars.maxthreads);
  

  // Test initialization
  blitz::Array<double, 4> ts_ls_exactdb_brdfunc_exp(tst_obj.ts_ls_exactdb_brdfunc().shape());
  ts_ls_exactdb_brdfunc_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_exactdb_brdfunc(), ts_ls_exactdb_brdfunc_exp, 1e-10);
  blitz::Array<double, 4> ts_ls_brdf_f_0_exp(tst_obj.ts_ls_brdf_f_0().shape());
  ts_ls_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_brdf_f_0(), ts_ls_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 4> ts_ls_brdf_f_exp(tst_obj.ts_ls_brdf_f().shape());
  ts_ls_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_brdf_f(), ts_ls_brdf_f_exp, 1e-10);
  blitz::Array<double, 4> ts_ls_user_brdf_f_0_exp(tst_obj.ts_ls_user_brdf_f_0().shape());
  ts_ls_user_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_user_brdf_f_0(), ts_ls_user_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 4> ts_ls_user_brdf_f_exp(tst_obj.ts_ls_user_brdf_f().shape());
  ts_ls_user_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_user_brdf_f(), ts_ls_user_brdf_f_exp, 1e-10);
  blitz::Array<double, 3> ts_ls_emissivity_exp(tst_obj.ts_ls_emissivity().shape());
  ts_ls_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_emissivity(), ts_ls_emissivity_exp, 1e-10);
  blitz::Array<double, 3> ts_ls_user_emissivity_exp(tst_obj.ts_ls_user_emissivity().shape());
  ts_ls_user_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_ls_user_emissivity(), ts_ls_user_emissivity_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_linsup_ss_atmos)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Linsup_Ss_Atmos tst_obj = Lidort_Linsup_Ss_Atmos();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_ss().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_ss().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_ss().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_ss().extent(3), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_db().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_db().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_columnwf_db().extent(2), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(3), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_ss().extent(4), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_db().extent(0), lid_pars.max_atmoswfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_db().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_db().extent(2), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_profilewf_db().extent(3), lid_pars.max_geometries);
  

  // Test initialization
  blitz::Array<double, 4> ts_columnwf_ss_exp(tst_obj.ts_columnwf_ss().shape());
  ts_columnwf_ss_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_columnwf_ss(), ts_columnwf_ss_exp, 1e-10);
  blitz::Array<double, 3> ts_columnwf_db_exp(tst_obj.ts_columnwf_db().shape());
  ts_columnwf_db_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_columnwf_db(), ts_columnwf_db_exp, 1e-10);
  blitz::Array<double, 5> ts_profilewf_ss_exp(tst_obj.ts_profilewf_ss().shape());
  ts_profilewf_ss_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_profilewf_ss(), ts_profilewf_ss_exp, 1e-10);
  blitz::Array<double, 4> ts_profilewf_db_exp(tst_obj.ts_profilewf_db().shape());
  ts_profilewf_db_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_profilewf_db(), ts_profilewf_db_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_linsup_ss_surf)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Linsup_Ss_Surf tst_obj = Lidort_Linsup_Ss_Surf();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf_db().extent(0), lid_pars.max_surfacewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf_db().extent(1), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_surfacewf_db().extent(2), lid_pars.max_geometries);
  

  // Test initialization
  blitz::Array<double, 3> ts_surfacewf_db_exp(tst_obj.ts_surfacewf_db().shape());
  ts_surfacewf_db_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_surfacewf_db(), ts_surfacewf_db_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_linsup_ss)
{
  // Test constructor
  Lidort_Linsup_Ss tst_obj = Lidort_Linsup_Ss();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(lidort_linsup_sleave)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Linsup_Sleave tst_obj = Lidort_Linsup_Sleave();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_isotropic().extent(0), lid_pars.max_sleavewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_isotropic().extent(1), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_userangles().extent(0), lid_pars.max_sleavewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_userangles().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_userangles().extent(2), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_userangles().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_f_0().extent(0), lid_pars.max_sleavewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_f_0().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_slterm_f_0().extent(3), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_user_slterm_f_0().extent(0), lid_pars.max_sleavewfs);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_user_slterm_f_0().extent(1), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_user_slterm_f_0().extent(2), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_lssl_user_slterm_f_0().extent(3), lid_pars.maxbeams);
  

  // Test initialization
  blitz::Array<double, 2> ts_lssl_slterm_isotropic_exp(tst_obj.ts_lssl_slterm_isotropic().shape());
  ts_lssl_slterm_isotropic_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lssl_slterm_isotropic(), ts_lssl_slterm_isotropic_exp, 1e-10);
  blitz::Array<double, 4> ts_lssl_slterm_userangles_exp(tst_obj.ts_lssl_slterm_userangles().shape());
  ts_lssl_slterm_userangles_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lssl_slterm_userangles(), ts_lssl_slterm_userangles_exp, 1e-10);
  blitz::Array<double, 4> ts_lssl_slterm_f_0_exp(tst_obj.ts_lssl_slterm_f_0().shape());
  ts_lssl_slterm_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lssl_slterm_f_0(), ts_lssl_slterm_f_0_exp, 1e-10);
  blitz::Array<double, 4> ts_lssl_user_slterm_f_0_exp(tst_obj.ts_lssl_user_slterm_f_0().shape());
  ts_lssl_user_slterm_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lssl_user_slterm_f_0(), ts_lssl_user_slterm_f_0_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_linsup_inout)
{
  // Test constructor
  Lidort_Linsup_Inout tst_obj = Lidort_Linsup_Inout();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(lidort_main_outputs)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Main_Outputs tst_obj = Lidort_Main_Outputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_intensity().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_intensity().extent(1), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_intensity().extent(2), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_intensity().extent(3), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_mean_intensity().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_mean_intensity().extent(1), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_mean_intensity().extent(2), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_mean_intensity().extent(3), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_integral().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_integral().extent(1), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_integral().extent(2), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_flux_integral().extent(3), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct().extent(1), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnflux_direct().extent(2), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmean_direct().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmean_direct().extent(1), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_dnmean_direct().extent(2), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_fourier_saved().extent(0), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_fourier_saved().extent(1), lid_pars.maxthreads);
  

  // Test initialization
  blitz::Array<double, 4> ts_intensity_exp(tst_obj.ts_intensity().shape());
  ts_intensity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_intensity(), ts_intensity_exp, 1e-10);
  blitz::Array<double, 4> ts_mean_intensity_exp(tst_obj.ts_mean_intensity().shape());
  ts_mean_intensity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_mean_intensity(), ts_mean_intensity_exp, 1e-10);
  blitz::Array<double, 4> ts_flux_integral_exp(tst_obj.ts_flux_integral().shape());
  ts_flux_integral_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_flux_integral(), ts_flux_integral_exp, 1e-10);
  blitz::Array<double, 3> ts_dnflux_direct_exp(tst_obj.ts_dnflux_direct().shape());
  ts_dnflux_direct_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_dnflux_direct(), ts_dnflux_direct_exp, 1e-10);
  blitz::Array<double, 3> ts_dnmean_direct_exp(tst_obj.ts_dnmean_direct().shape());
  ts_dnmean_direct_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_dnmean_direct(), ts_dnmean_direct_exp, 1e-10);
  blitz::Array<int, 2> ts_fourier_saved_exp(tst_obj.ts_fourier_saved().shape());
  ts_fourier_saved_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_fourier_saved(), ts_fourier_saved_exp, 1e-10);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_geometries(), 0);
  
}

BOOST_AUTO_TEST_CASE(lidort_exception_handling)
{
  // Test constructor
  Lidort_Exception_Handling tst_obj = Lidort_Exception_Handling();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_status_inputcheck(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_ncheckmessages(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_status_calculation(), 0);
  
}

BOOST_AUTO_TEST_CASE(lidort_input_exception_handling)
{
  // Test constructor
  Lidort_Input_Exception_Handling tst_obj = Lidort_Input_Exception_Handling();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_status_inputread(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_ninputmessages(), 0);
  
}

BOOST_AUTO_TEST_CASE(lidort_outputs)
{
  // Test constructor
  Lidort_Outputs tst_obj = Lidort_Outputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(lidort_sup_brdf)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Sup_Brdf tst_obj = Lidort_Sup_Brdf();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_exactdb_brdfunc().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_exactdb_brdfunc().extent(1), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_exactdb_brdfunc().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f_0().extent(1), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f_0().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f().extent(1), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_brdf_f().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f_0().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f_0().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_brdf_f().extent(2), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_emissivity().extent(0), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_emissivity().extent(1), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_emissivity().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_emissivity().extent(1), lid_pars.maxthreads);
  

  // Test initialization
  blitz::Array<double, 3> ts_exactdb_brdfunc_exp(tst_obj.ts_exactdb_brdfunc().shape());
  ts_exactdb_brdfunc_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_exactdb_brdfunc(), ts_exactdb_brdfunc_exp, 1e-10);
  blitz::Array<double, 3> ts_brdf_f_0_exp(tst_obj.ts_brdf_f_0().shape());
  ts_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_brdf_f_0(), ts_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 3> ts_brdf_f_exp(tst_obj.ts_brdf_f().shape());
  ts_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_brdf_f(), ts_brdf_f_exp, 1e-10);
  blitz::Array<double, 3> ts_user_brdf_f_0_exp(tst_obj.ts_user_brdf_f_0().shape());
  ts_user_brdf_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_brdf_f_0(), ts_user_brdf_f_0_exp, 1e-10);
  blitz::Array<double, 3> ts_user_brdf_f_exp(tst_obj.ts_user_brdf_f().shape());
  ts_user_brdf_f_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_brdf_f(), ts_user_brdf_f_exp, 1e-10);
  blitz::Array<double, 2> ts_emissivity_exp(tst_obj.ts_emissivity().shape());
  ts_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_emissivity(), ts_emissivity_exp, 1e-10);
  blitz::Array<double, 2> ts_user_emissivity_exp(tst_obj.ts_user_emissivity().shape());
  ts_user_emissivity_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_emissivity(), ts_user_emissivity_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_sup_sleave)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Sup_Sleave tst_obj = Lidort_Sup_Sleave();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_isotropic().extent(0), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_userangles().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_userangles().extent(1), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_userangles().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_f_0().extent(1), lid_pars.maxstreams);
  BOOST_CHECK_EQUAL(tst_obj.ts_slterm_f_0().extent(2), lid_pars.maxbeams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_slterm_f_0().extent(0), lid_pars.maxmoments+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_slterm_f_0().extent(1), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_slterm_f_0().extent(2), lid_pars.maxbeams);
  

  // Test initialization
  blitz::Array<double, 1> ts_slterm_isotropic_exp(tst_obj.ts_slterm_isotropic().shape());
  ts_slterm_isotropic_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_slterm_isotropic(), ts_slterm_isotropic_exp, 1e-10);
  blitz::Array<double, 3> ts_slterm_userangles_exp(tst_obj.ts_slterm_userangles().shape());
  ts_slterm_userangles_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_slterm_userangles(), ts_slterm_userangles_exp, 1e-10);
  blitz::Array<double, 3> ts_slterm_f_0_exp(tst_obj.ts_slterm_f_0().shape());
  ts_slterm_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_slterm_f_0(), ts_slterm_f_0_exp, 1e-10);
  blitz::Array<double, 3> ts_user_slterm_f_0_exp(tst_obj.ts_user_slterm_f_0().shape());
  ts_user_slterm_f_0_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_slterm_f_0(), ts_user_slterm_f_0_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_sup_ss)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Sup_Ss tst_obj = Lidort_Sup_Ss();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_intensity_ss().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_intensity_ss().extent(1), lid_pars.max_geometries);
  BOOST_CHECK_EQUAL(tst_obj.ts_intensity_ss().extent(2), lid_pars.max_directions);
  BOOST_CHECK_EQUAL(tst_obj.ts_intensity_db().extent(0), lid_pars.max_user_levels);
  BOOST_CHECK_EQUAL(tst_obj.ts_intensity_db().extent(1), lid_pars.max_geometries);
  

  // Test initialization
  blitz::Array<double, 3> ts_intensity_ss_exp(tst_obj.ts_intensity_ss().shape());
  ts_intensity_ss_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_intensity_ss(), ts_intensity_ss_exp, 1e-10);
  blitz::Array<double, 2> ts_intensity_db_exp(tst_obj.ts_intensity_db().shape());
  ts_intensity_db_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_intensity_db(), ts_intensity_db_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_sup_inout)
{
  // Test constructor
  Lidort_Sup_Inout tst_obj = Lidort_Sup_Inout();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(lidort_fixed_boolean)
{
  // Test constructor
  Lidort_Fixed_Boolean tst_obj = Lidort_Fixed_Boolean();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_do_fullrad_mode(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_sscorr_truncation(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_ss_external(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_ssfull(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_thermal_emission(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_surface_emission(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_plane_parallel(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_brdf_surface(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_upwelling(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_dnwelling(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_surface_leaving(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_sl_isotropic(), false);
  
}

BOOST_AUTO_TEST_CASE(lidort_fixed_control)
{
  // Test constructor
  Lidort_Fixed_Control tst_obj = Lidort_Fixed_Control();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_nstreams(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_nlayers(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_nfinelayers(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_thermal_coeffs(), 0);
  BOOST_CHECK_CLOSE(tst_obj.ts_lidort_accuracy(), 0, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_fixed_sunrays)
{
  // Test constructor
  Lidort_Fixed_Sunrays tst_obj = Lidort_Fixed_Sunrays();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_CLOSE(tst_obj.ts_flux_factor(), 0, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_fixed_uservalues)
{
  // Test constructor
  Lidort_Fixed_Uservalues tst_obj = Lidort_Fixed_Uservalues();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_n_user_streams(), 0);
  BOOST_CHECK_EQUAL(tst_obj.ts_n_user_levels(), 0);
  
}

BOOST_AUTO_TEST_CASE(lidort_fixed_chapman)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Fixed_Chapman tst_obj = Lidort_Fixed_Chapman();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_height_grid().extent(0), lid_pars.maxlayers+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_pressure_grid().extent(0), lid_pars.maxlayers+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_temperature_grid().extent(0), lid_pars.maxlayers+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_finegrid().extent(0), lid_pars.maxlayers);
  

  // Test initialization
  blitz::Array<double, 1> ts_height_grid_exp(tst_obj.ts_height_grid().shape());
  ts_height_grid_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_height_grid(), ts_height_grid_exp, 1e-10);
  blitz::Array<double, 1> ts_pressure_grid_exp(tst_obj.ts_pressure_grid().shape());
  ts_pressure_grid_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_pressure_grid(), ts_pressure_grid_exp, 1e-10);
  blitz::Array<double, 1> ts_temperature_grid_exp(tst_obj.ts_temperature_grid().shape());
  ts_temperature_grid_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_temperature_grid(), ts_temperature_grid_exp, 1e-10);
  blitz::Array<int, 1> ts_finegrid_exp(tst_obj.ts_finegrid().shape());
  ts_finegrid_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_finegrid(), ts_finegrid_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_rfindex_parameter(), 0, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_fixed_optical)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Fixed_Optical tst_obj = Lidort_Fixed_Optical();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_deltau_vert_input().extent(0), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_deltau_vert_input().extent(1), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_phasmoms_total_input().extent(0), lid_pars.maxmoments_input+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_phasmoms_total_input().extent(1), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_phasmoms_total_input().extent(2), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_thermal_bb_input().extent(0), lid_pars.maxlayers+1);
  BOOST_CHECK_EQUAL(tst_obj.ts_thermal_bb_input().extent(1), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_lambertian_albedo().extent(0), lid_pars.maxthreads);
  BOOST_CHECK_EQUAL(tst_obj.ts_surface_bb_input().extent(0), lid_pars.maxthreads);
  

  // Test initialization
  blitz::Array<double, 2> ts_deltau_vert_input_exp(tst_obj.ts_deltau_vert_input().shape());
  ts_deltau_vert_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_deltau_vert_input(), ts_deltau_vert_input_exp, 1e-10);
  blitz::Array<double, 3> ts_phasmoms_total_input_exp(tst_obj.ts_phasmoms_total_input().shape());
  ts_phasmoms_total_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_phasmoms_total_input(), ts_phasmoms_total_input_exp, 1e-10);
  blitz::Array<double, 2> ts_thermal_bb_input_exp(tst_obj.ts_thermal_bb_input().shape());
  ts_thermal_bb_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_thermal_bb_input(), ts_thermal_bb_input_exp, 1e-10);
  blitz::Array<double, 1> ts_lambertian_albedo_exp(tst_obj.ts_lambertian_albedo().shape());
  ts_lambertian_albedo_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_lambertian_albedo(), ts_lambertian_albedo_exp, 1e-10);
  blitz::Array<double, 1> ts_surface_bb_input_exp(tst_obj.ts_surface_bb_input().shape());
  ts_surface_bb_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_surface_bb_input(), ts_surface_bb_input_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_fixed_inputs)
{
  // Test constructor
  Lidort_Fixed_Inputs tst_obj = Lidort_Fixed_Inputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}

BOOST_AUTO_TEST_CASE(lidort_modified_boolean)
{
  // Test constructor
  Lidort_Modified_Boolean tst_obj = Lidort_Modified_Boolean();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_do_sscorr_nadir(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_sscorr_outgoing(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_double_convtest(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_solar_sources(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_refractive_geometry(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_chapman_function(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_rayleigh_only(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_isotropic_only(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_no_azimuth(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_all_fourier(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_deltam_scaling(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_solution_saving(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_bvp_telescoping(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_user_streams(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_additional_mvout(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_mvout_only(), false);
  BOOST_CHECK_EQUAL(tst_obj.ts_do_thermal_transonly(), false);
  
}

BOOST_AUTO_TEST_CASE(lidort_modified_control)
{
  // Test constructor
  Lidort_Modified_Control tst_obj = Lidort_Modified_Control();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_nmoments_input(), 0);
  
}

BOOST_AUTO_TEST_CASE(lidort_modified_sunrays)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Modified_Sunrays tst_obj = Lidort_Modified_Sunrays();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_beam_szas().extent(0), lid_pars.maxbeams);
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_nbeams(), 0);
  blitz::Array<double, 1> ts_beam_szas_exp(tst_obj.ts_beam_szas().shape());
  ts_beam_szas_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_beam_szas(), ts_beam_szas_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_modified_uservalues)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Modified_Uservalues tst_obj = Lidort_Modified_Uservalues();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_user_relazms().extent(0), lid_pars.max_user_relazms);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_angles_input().extent(0), lid_pars.max_user_streams);
  BOOST_CHECK_EQUAL(tst_obj.ts_user_levels().extent(0), lid_pars.max_user_levels);
  

  // Test initialization
  BOOST_CHECK_EQUAL(tst_obj.ts_n_user_relazms(), 0);
  blitz::Array<double, 1> ts_user_relazms_exp(tst_obj.ts_user_relazms().shape());
  ts_user_relazms_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_relazms(), ts_user_relazms_exp, 1e-10);
  blitz::Array<double, 1> ts_user_angles_input_exp(tst_obj.ts_user_angles_input().shape());
  ts_user_angles_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_angles_input(), ts_user_angles_input_exp, 1e-10);
  blitz::Array<double, 1> ts_user_levels_exp(tst_obj.ts_user_levels().shape());
  ts_user_levels_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_user_levels(), ts_user_levels_exp, 1e-10);
  BOOST_CHECK_CLOSE(tst_obj.ts_geometry_specheight(), 0, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_modified_chapman)
{
  // Test constructor
  Lidort_Modified_Chapman tst_obj = Lidort_Modified_Chapman();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  BOOST_CHECK_CLOSE(tst_obj.ts_earth_radius(), 0, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_modified_optical)
{
  // Used for checking dimensions
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Test constructor
  Lidort_Modified_Optical tst_obj = Lidort_Modified_Optical();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  BOOST_CHECK_EQUAL(tst_obj.ts_omega_total_input().extent(0), lid_pars.maxlayers);
  BOOST_CHECK_EQUAL(tst_obj.ts_omega_total_input().extent(1), lid_pars.maxthreads);
  

  // Test initialization
  blitz::Array<double, 2> ts_omega_total_input_exp(tst_obj.ts_omega_total_input().shape());
  ts_omega_total_input_exp = 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.ts_omega_total_input(), ts_omega_total_input_exp, 1e-10);
  
}

BOOST_AUTO_TEST_CASE(lidort_modified_inputs)
{
  // Test constructor
  Lidort_Modified_Inputs tst_obj = Lidort_Modified_Inputs();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  

  // Test initialization
  
}



BOOST_AUTO_TEST_SUITE_END()
       
