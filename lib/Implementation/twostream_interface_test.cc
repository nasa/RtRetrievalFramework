#include "twostream_interface.h"
#include "unit_test_support.h"
#include "lidort_interface_types.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(twostream_interface, GlobalFixture)

BOOST_AUTO_TEST_CASE(twostream_l_master)
{

  // For obtaining maximum sizes
  Lidort_Pars lid_pars = Lidort_Pars::instance();
 
  // Dimensioning integers
  int nbeams         = 2;
  int n_user_angles  = 2;
  int n_user_relazms = 2;
  int nlayers        = 3;
  int npars          = 3;
  int nspars         = 1;
  int nstreams_brdf  = 50;
  int nbrdf          = 4;

  int n_geometries   = nbeams * n_user_angles*n_user_relazms;
  int ntotal         = 2*nlayers;

  // Directional Flags
  bool do_upwelling = true;
  bool do_dnwelling = true;

  // Sources control, including thermal   
  bool do_thermal_emission;
  bool do_surface_emission;
  bool do_solar_sources;
    
  // Plane parallel and deltam-2stream scaling flags
  bool do_plane_parallel, do_delta2s_scaling;

  // BRDF surface flag
  bool do_brdf_surface;

  // Linearization flags
  bool do_profile_wfs, do_surface_wfs;

  // Linearization control
  Array<bool, 1> layer_vflag( nlayers );
  Array<int, 1> layer_vnumber( nlayers );

  int n_surface_wfs;

  // Geometry
  Array<double, 1> beam_szas( nbeams );
  Array<double, 1> user_angles( n_user_angles );
  Array<double, 1> user_relazms( n_user_relazms );

  // Stream value
  double stream_value;

  // Cox-Munk Surface control (threaded)
  double wind_speed;

  // BRDF Fourier components (NOT threaded)
  // 0 and 1 Fourier components of BRDF, following order (same all threads)
  //   incident solar directions,  reflected quadrature stream
  //   incident quadrature stream, reflected quadrature stream
  //   incident solar directions,  reflected user streams
  //   incident quadrature stream, reflected user streams
  Array<double, 2> brdf_f_0(2, nbeams);  // ( 0:1, NBEAMS )

  Array<double, 1> brdf_f(2); // ( 0:1 )
  Array<double, 2> ubrdf_f(2, n_user_angles); // ( 0:1, N_USER_ANGLES )

  // Linearized BRDF fourier components
  // 0 and 1 Fourier components of BRDF, following order (same all threads)
  //   incident solar directions,  reflected quadrature stream
  //   incident quadrature stream, reflected quadrature stream
  //   incident solar directions,  reflected user streams
  //   incident quadrature stream, reflected user streams
  Array<double, 3> ls_brdf_f_0(nspars, 2, nbeams); // (NSPARS,0:1,NBEAMS)
  Array<double, 2> ls_brdf_f(nspars, 2); // (NSPARS,0:1)
  Array<double, 3> ls_ubrdf_f(nspars, 2, n_user_angles); // (NSPARS,0:1,N_USER_ANGLES)

  // Other BRDF variables
  bool do_shadow_effect;
  Array<bool, 1> lambertian_kernel_flag(nbrdf);
  Array<bool, 1> do_kernel_factor_wfs(nbrdf);
  Array<bool, 2> do_kernel_params_wfs(nbrdf,nbrdf);
  Array<bool, 1> do_kparams_derivs(nbrdf);
  int n_brdf_kernels;
  Array<int, 1> n_brdf_parameters(nbrdf);
  Array<int, 1> which_brdf(nbrdf);
  Array<double, 1> brdf_factors(nbrdf);
  Array<double, 2> brdf_parameters(nbrdf,nbrdf);

  // Thermal variables
  double surfbb = 0;
  Array<double, 1> thermal_bb_input(nlayers+1); // ( 0:NLAYERS )
  double emissivity;
  Array<double, 1> ls_emissivity(nspars);

  // Flux factor
  double flux_factor;

  // Height and earth radius
  double earth_radius;
  Array<double, 1> height_grid(nlayers+1); // ( 0:NLAYERS )

  // Atmospheric Optical properties
  Array<double, 1> deltau_input(nlayers);
  Array<double, 1> omega_input(nlayers);
  Array<double, 1> asymm_input(nlayers);
  Array<double, 1> d2s_scaling(nlayers);

  // Linearized optical properties
  Array<double, 2> l_deltau_input(nlayers, npars);
  Array<double, 2> l_omega_input(nlayers, npars);
  Array<double, 2> l_asymm_input(nlayers, npars);
  Array<double, 2> l_d2s_scaling(nlayers, npars);

  // Results
  Array<double, 1> intensity_toa(n_geometries);
  Array<double, 1> intensity_boa(n_geometries);
  
  Array<double, 3> profilewf_toa(n_geometries, nlayers, npars);
  Array<double, 3> profilewf_boa(n_geometries, nlayers, npars);

  Array<double, 2> surfacewf_toa(n_geometries, nspars);
  Array<double, 2> surfacewf_boa(n_geometries, nspars);

  // Exception handling
  int status_inputcheck;
  int c_nmessages;
  std::vector<std::string> c_messages;
  std::vector<std::string> c_actions;

  int status_execution;
  std::string e_message, e_trace_1, e_trace_2;

  //  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //          OTHER ARGUMENTS
  //  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  //  Other flags
  bool do_fullquadrature;

  // Other variables
  Array<double, 1> abc0(3), sca0(3), asy0(3);
  Array<double, 1> abc(3), sca(3);
  double column;
  double omega, omega1, s2, m1, m2;

  Array<double, 1> abc_save(3);
  abc_save = 0.05, 0.02, 0.01;
  Array<double, 1> sca_save(3);
  sca_save = 0.15, 0.28, 0.34;
  Array<double, 1> asy_save(3);
  asy_save = 0.70, 0.76, 0.77;

  Array<double, 1> abc2_save(3);
  abc2_save = 0.04, 0.03, 0.03;
  Array<double, 1> sca2_save(3);
  sca2_save = 0.14, 0.35, 0.41;
  Array<double, 1> b2(3);
  b2 = 0.5, 0.5, 0.5;

  double column_save = 1.0e0;
  wind_speed = 10.e0;

  // Set up inputs
  do_brdf_surface     = false;
  do_fullquadrature   = false;
  do_delta2s_scaling  = true;
  do_solar_sources    = true;
  do_thermal_emission = false;
  do_surface_emission = false;
  do_plane_parallel   = false;

  // Set DO_BRDF_SURFACE to true
  do_brdf_surface = true;

  //  Other 2-stream flags

  //  Flux factor
  flux_factor = 1.0e0;

  beam_szas = 57.0e0, 66.0e0;
  user_angles = 35.0e0, 55.0e0;
  user_relazms= 10.0e0, 170.0e0;

  earth_radius = 6371.0e0;

  height_grid = 36.0e0, 20.0e0, 10.0e0, 0.0e0;

  // Two choices of stream value................ CHOOSE One !!!!!
  if ( do_fullquadrature )
    stream_value = sqrt(1.0e0 / 3.0e0 );
  else
    stream_value = 0.5e0;

  // Set up the BRDF stuff, Cox Munk
  // ======================================================

  do_solar_sources          = true;
  do_shadow_effect          = true;
  lambertian_kernel_flag(0) = false;
  do_kernel_factor_wfs      = false;
  do_kernel_params_wfs      = false;
  do_kernel_params_wfs(0,0) = true;
  do_kparams_derivs         = false;
  n_brdf_kernels            = 1;
  n_brdf_parameters         = 0;
  n_brdf_parameters(0)      = 3;
  which_brdf                = 0;
  which_brdf(0)             = 9;
  n_surface_wfs             = 0;
  brdf_factors              = 0.e0;
  brdf_factors(0)           = 1.e0;
  brdf_parameters           = 0.e0;
  brdf_parameters(0,1)      = 1.33e0 * 1.33e0; // Square of refractive index
  brdf_parameters(0,0)      = 0.003e0 + 0.00512e0 * wind_speed; // Slope squared

  // Call linearized BRDF supplement
  Twostream_Ls_Brdf_Supplement twostream_brdf(
    // Max sizes
    nbeams, n_user_angles, lid_pars.max_user_obsgeoms, nstreams_brdf,
    lid_pars.max_brdf_kernels, lid_pars.max_brdf_parameters, nspars,
    // Instance sizes
    nbeams, n_user_angles, nstreams_brdf);

  // Inputs

  twostream_brdf.do_solar_sources(do_solar_sources);
  twostream_brdf.do_shadow_effect(do_shadow_effect);
  twostream_brdf.do_surface_emission(do_surface_emission);
  twostream_brdf.n_brdf_kernels(n_brdf_kernels);
  twostream_brdf.n_brdf_parameters(n_brdf_parameters);
       
  twostream_brdf.which_brdf(which_brdf);
  twostream_brdf.lambertian_kernel_flag(lambertian_kernel_flag);
                  
  twostream_brdf.brdf_factors(brdf_factors);
  twostream_brdf.brdf_parameters(brdf_parameters);
                          
  twostream_brdf.do_kernel_factor_wfs(do_kernel_factor_wfs);
  twostream_brdf.do_kernel_params_wfs(do_kernel_params_wfs);

  twostream_brdf.beam_szas(beam_szas);
  twostream_brdf.user_angles(user_angles);
  //twostream_brdf.user_relazms(user_relazms);
  twostream_brdf.stream_value(stream_value);

  // Make calculation
  twostream_brdf.run();

  // Outputs                  
  do_kparams_derivs = twostream_brdf.do_kparams_derivs();
  n_surface_wfs = twostream_brdf.n_surface_wfs();
  twostream_brdf.n_kernel_factor_wfs();
  twostream_brdf.n_kernel_params_wfs();

 brdf_f_0.reference(twostream_brdf.brdf_f_0());
  brdf_f.reference(twostream_brdf.brdf_f());
  ubrdf_f.reference(twostream_brdf.ubrdf_f());

  ls_brdf_f_0.reference(twostream_brdf.ls_brdf_f_0());
  ls_brdf_f.reference(twostream_brdf.ls_brdf_f());
  ls_ubrdf_f.reference(twostream_brdf.ls_ubrdf_f());

  emissivity = twostream_brdf.emissivity();
  ls_emissivity.reference(twostream_brdf.ls_emissivity());

  // Baseline calculation 2 : RADIANCE+ PROFILE/SURFACE WFS
  // ======================================================

  //  Control for 3 profile WFs
  do_profile_wfs = true;
  do_surface_wfs = true;
  //do_sim_only    = false;

  n_surface_wfs = 1;
  for(int n = 0; n < nlayers; n++) {
    layer_vflag(n)   = true;
    layer_vnumber(n) = 3;
  }

  //  Optical property inputs
  column = column_save;
  for(int n = 0; n < nlayers; n++) {
    abc(n) = column * abc_save(n) + abc2_save(n);
    sca(n) = sca_save(n) + sca2_save(n);
    deltau_input(n) = abc(n) + sca(n);
    omega_input(n)  = sca(n) /  deltau_input(n);
    asymm_input(n)  = asy_save(n)* sca2_save(n) / sca(n);
    m1 = b2(n) * sca_save(n);
    m2 = 5.0e0 * asy_save(n) * asy_save(n) * sca2_save(n);
    d2s_scaling(n)  = ( m1 + m2 ) / sca(n) / 5.0e0;
  }

  //  First WF is w.r.t abc_save(n)
  for(int n = 0; n < nlayers; n++) {
    omega = omega_input(n);
    l_deltau_input(n,0) = column;
    l_omega_input (n,0) = - column * omega /deltau_input(n);
    l_asymm_input (n,0) = 0.0e0;
    l_d2s_scaling (n,0) = 0.0e0;
  }

  //  Second WF is w.r.t sca_save(n)
  for(int n = 0; n < nlayers; n++) {
    sca(n) = sca_save(n) + sca2_save(n);
    s2     =  sca2_save(n) / sca(n);
    m1     = b2(n) * sca_save(n) / sca(n);
    m2     = 5.0e0 * asy_save(n) * asy_save(n) * s2;
    omega1 = 1.0e0-omega_input(n);
    l_deltau_input(n,1) = 1.0e0;
    l_omega_input (n,1) = omega1 / deltau_input(n);
    l_asymm_input (n,1) = - asymm_input(n) / sca(n);
    l_d2s_scaling (n,1) = ( b2(n) - (m1+m2) ) / sca(n) / 5.0e0;
  }

  //  Third WF is w.r.t. asy_save(n)
  for(int n = 0; n < nlayers; n++) {
    sca(n) = sca_save(n) + sca2_save(n);
    s2     =  sca2_save(n) / sca(n);
    l_deltau_input(n,2) = 0.0e0;
    l_omega_input (n,2) = 0.0e0;
    l_asymm_input (n,2) = s2;
    l_d2s_scaling (n,2) = 2.0e0 * asy_save(n) * s2;
  }

  //  Normalize, profile WFs
  for(int n = 0; n < nlayers; n++) {
    l_deltau_input(n,0) = l_deltau_input(n,0) * abc_save(n);
    l_omega_input (n,0) = l_omega_input (n,0) * abc_save(n);
    l_asymm_input (n,0) = l_asymm_input (n,0) * abc_save(n);
    l_d2s_scaling (n,0) = l_d2s_scaling (n,0) * abc_save(n);
    l_deltau_input(n,1) = l_deltau_input(n,1) * sca_save(n);
    l_omega_input (n,1) = l_omega_input (n,1) * sca_save(n);
    l_asymm_input (n,1) = l_asymm_input (n,1) * sca_save(n);
    l_d2s_scaling (n,1) = l_d2s_scaling (n,1) * sca_save(n);
    l_deltau_input(n,2) = l_deltau_input(n,2) * asy_save(n);
    l_omega_input (n,2) = l_omega_input (n,2) * asy_save(n);
    l_asymm_input (n,2) = l_asymm_input (n,2) * asy_save(n);
    l_d2s_scaling (n,2) = l_d2s_scaling (n,2) * asy_save(n);
  }

  Twostream_Lps_Master twostream_rt = Twostream_Lps_Master( 
    // Max sizes
    nlayers, lid_pars.maxtotal, lid_pars.max_messages,
    nbeams, n_geometries, 
    n_user_angles, n_user_relazms, lid_pars.max_user_obsgeoms, 
    npars, nspars, lid_pars.max_sleavewfs, 
    // Instance sizes
    nlayers, ntotal, n_user_angles, n_user_relazms, nbeams, earth_radius, n_geometries);
      
  // Input
  twostream_rt.do_upwelling(do_upwelling);
  twostream_rt.do_dnwelling(do_dnwelling);
  twostream_rt.do_plane_parallel(do_plane_parallel);

  twostream_rt.do_d2s_scaling(do_delta2s_scaling);
  twostream_rt.do_brdf_surface(do_brdf_surface);

  twostream_rt.do_solar_sources(do_solar_sources);
  twostream_rt.do_thermal_emission(do_thermal_emission);
  twostream_rt.do_surface_emission(do_surface_emission);

  twostream_rt.do_profile_wfs(do_profile_wfs);
  twostream_rt.do_surface_wfs(do_surface_wfs);
  //twostream_rt.do_sim_only(do_sim_only);

  twostream_rt.layer_vary_flag(layer_vflag);
  twostream_rt.layer_vary_number(layer_vnumber);

  twostream_rt.n_surface_wfs(n_surface_wfs);

  twostream_rt.beam_szas(beam_szas);
  twostream_rt.user_angles(user_angles);
  twostream_rt.user_relazms(user_relazms);
  twostream_rt.stream_value(stream_value);

  twostream_rt.bvpscalefactor(1.0);

  twostream_rt.brdf_f_0(brdf_f_0);
  twostream_rt.brdf_f(brdf_f);
  twostream_rt.ubrdf_f(ubrdf_f);

  twostream_rt.thermal_bb_input(thermal_bb_input);
  twostream_rt.surfbb(surfbb);
  twostream_rt.emissivity(emissivity);
  twostream_rt.ls_emissivity(ls_emissivity);

  twostream_rt.ls_brdf_f_0(ls_brdf_f_0);
  twostream_rt.ls_brdf_f(ls_brdf_f);
  twostream_rt.ls_ubrdf_f(ls_ubrdf_f);

  twostream_rt.flux_factor(flux_factor);

  twostream_rt.height_grid(height_grid);
  twostream_rt.deltau_input(deltau_input);
  twostream_rt.omega_input(omega_input);

  twostream_rt.asymm_input(asymm_input);
  twostream_rt.d2s_scaling(d2s_scaling);

  twostream_rt.l_deltau_input(l_deltau_input);
  twostream_rt.l_omega_input(l_omega_input);
  twostream_rt.l_asymm_input(l_asymm_input);
  twostream_rt.l_d2s_scaling(l_d2s_scaling);

  // Run radiative transfer
  twostream_rt.run();

  // Output
  intensity_toa.reference(twostream_rt.intensity_toa());
  profilewf_toa.reference(twostream_rt.profilewf_toa());
  surfacewf_toa.reference(twostream_rt.surfacewf_toa());
  
  intensity_boa.reference(twostream_rt.intensity_boa());
  profilewf_boa.reference(twostream_rt.profilewf_boa());
  surfacewf_boa.reference(twostream_rt.surfacewf_boa());
  
  status_inputcheck = twostream_rt.status_inputcheck();
  c_nmessages = twostream_rt.c_nmessages();
  c_messages = twostream_rt.c_messages();
  c_actions = twostream_rt.c_actions();
  status_execution = twostream_rt.status_execution();
  e_message = twostream_rt.e_message();
  e_trace_1 = twostream_rt.e_trace_1();
  e_trace_2 = twostream_rt.e_trace_2();

  //  Exception handling
  if (status_inputcheck == 1) {
    std::cerr << "INPUT Check failed from Baseline Run # 1" << std::endl
	      << " - Number of Messages = " << c_nmessages << std::endl;
    for(int k = 0; k < c_nmessages; k++) {
      std::cerr << " - Message # " << k << ": " << c_messages[k] << std::endl
		<< " - Action  # " << k << ": " << c_actions[k] << std::endl;
    }
  }
  if (status_execution == 1) {
    std::cerr << "EXECUTION failed from Baseline Run # 1" << std::endl
	      << " - Print 1 Message and 2 Traces" << std::endl
	      << e_message << std::endl
	      << e_trace_1 << std::endl
	      << e_trace_2 << std::endl;
  }

  Array<double, 2> intensity_toa_expt;
  Array<double, 2> intensity_boa_expt;
  Array<double, 4> profilewf_toa_expt;
  Array<double, 4> profilewf_boa_expt;
  Array<double, 3> surfacewf_toa_expt;
  Array<double, 3> surfacewf_boa_expt;

  IfstreamCs twostream_expected(test_data_dir() + "expected/twostream/twostream_l_master");
  twostream_expected >> intensity_toa_expt
		     >> intensity_boa_expt
		     >> profilewf_toa_expt
		     >> profilewf_boa_expt
		     >> surfacewf_toa_expt
		     >> surfacewf_boa_expt;

  // Ordering of arrays are different so we need to compare one-dimensionally
  Range all = Range::all();
  BOOST_CHECK_MATRIX_CLOSE(intensity_toa(all), intensity_toa_expt(all, 0));
  BOOST_CHECK_MATRIX_CLOSE(intensity_boa(all), intensity_boa_expt(all, 0));
  for(int g = 0; g < n_geometries; g++) {
    for(int n = 0; n < nlayers; n++) {
      BOOST_CHECK_MATRIX_CLOSE(profilewf_toa(g, n, all), profilewf_toa_expt(g, n, all, 0));
      BOOST_CHECK_MATRIX_CLOSE(profilewf_boa(g, n, all), profilewf_boa_expt(g, n, all, 0));
    }
    BOOST_CHECK_MATRIX_CLOSE(surfacewf_boa(g, all), surfacewf_boa_expt(g, all, 0));
    BOOST_CHECK_MATRIX_CLOSE(surfacewf_toa(g, all), surfacewf_toa_expt(g, all, 0));
  }

  // Display result
  if(false) {
    std::cerr << std::scientific << std::setprecision(10);
    for(int g = 0; g < n_geometries; g++) {
      std::cerr << "intensity_toa: " << g << " = " << intensity_toa(g) << std::endl
		<< "intensity_boa: " << g << " = " << intensity_boa(g) << std::endl;
      for(int n = 0; n < nlayers; n++) {
	std::cerr << "profilewf_toa: " << g << " " << n << " = " << profilewf_toa(g, n, all, 0) << std::endl
		  << "profilewf_boa: " << g << " " << n << " = " << profilewf_boa(g, n, all, 0) << std::endl;
      }
    }
  }

}

BOOST_AUTO_TEST_SUITE_END()
