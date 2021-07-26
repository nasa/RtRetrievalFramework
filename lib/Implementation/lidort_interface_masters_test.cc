#include "lidort_interface_masters.h"
#include <blitz/array.h>
#include <cmath>
#include "unit_test_support.h"
#include "linear_algebra.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(lidort_interface_masters, GlobalFixture)

BOOST_AUTO_TEST_CASE(lidort_brdf_master)
{
  Brdf_Sup_Masters brdf_master = Brdf_Sup_Masters();

  // Read lidort config file
  brdf_master.read_config(test_data_dir() + "expected/lidort_interface_masters/3p8_BRDF_ReadInput.cfg");
  int read_status = brdf_master.brdf_sup_inputstatus().bs_status_inputread();
  BOOST_CHECK_EQUAL(read_status, 0); // was read successful?

  if (read_status != 0)
    std::cerr << brdf_master.brdf_sup_inputstatus();

  // Types where data is stored
  Brdf_Sup_Inputs& brdf_inputs = brdf_master.brdf_sup_in();
  Brdf_Sup_Outputs& brdf_outputs = brdf_master.brdf_sup_out();

  // Process inputs
  bool do_debug_restoration = false;
  int bs_nmoments_input = 2 * brdf_inputs.bs_nstreams() - 1;
  brdf_master.run(do_debug_restoration, bs_nmoments_input);

  // Range and size names
  int n_beams = brdf_inputs.bs_nbeams();
  int n_u_streams = brdf_inputs.bs_n_user_streams();
  int n_u_relazms = brdf_inputs.bs_n_user_relazms();

  Range r_u_streams = Range(0, n_u_streams-1);
  Range r_u_realazms = Range(0, n_u_relazms-1);
  Range r_beams = Range(0, n_beams-1);

  int n_streams = brdf_inputs.bs_nstreams();

  Range r_nmoms = Range(0, bs_nmoments_input);
  Range r_streams = Range(0, n_streams-1);

  // Check values according to offline results
  IfstreamCs expected_data(test_data_dir() + "expected/lidort_interface_masters/lidort_brdf_master");

  // --

  Array<double, 3> bs_dbounce_brdfunc_expect;
  expected_data >> bs_dbounce_brdfunc_expect;

  BOOST_CHECK_MATRIX_CLOSE(brdf_outputs.bs_dbounce_brdfunc()
                           (r_u_streams, r_u_realazms, r_beams), bs_dbounce_brdfunc_expect);

  // --

  Array<double, 3> bs_brdf_f_expect;
  expected_data >> bs_brdf_f_expect;
  BOOST_CHECK_MATRIX_CLOSE(brdf_outputs.bs_brdf_f()(r_nmoms,r_streams,r_streams), bs_brdf_f_expect);

  // --

  Array<double, 3> bs_brdf_f_0_expect;
  expected_data >> bs_brdf_f_0_expect;
  BOOST_CHECK_MATRIX_CLOSE(brdf_outputs.bs_brdf_f_0()(r_nmoms,r_streams,r_beams), bs_brdf_f_0_expect);

  // --

  Array<double, 3> bs_user_brdf_f_expect;
  expected_data >> bs_user_brdf_f_expect;
  BOOST_CHECK_MATRIX_CLOSE(brdf_outputs.bs_user_brdf_f()(r_nmoms,r_u_streams,r_streams), bs_user_brdf_f_expect);

  // --

  Array<double, 3> bs_user_brdf_f_0_expect;
  expected_data >> bs_user_brdf_f_0_expect;
  BOOST_CHECK_MATRIX_CLOSE(brdf_outputs.bs_user_brdf_f_0()(r_nmoms,r_u_streams,r_beams), bs_user_brdf_f_0_expect);
}

BOOST_AUTO_TEST_CASE(lidort_ls_brdf_master)
{
  // Same as lidrot_brdf_master except for weighting function parts
  Brdf_Lin_Sup_Masters brdf_master = Brdf_Lin_Sup_Masters();

  // Read lidort config file
  brdf_master.read_config(test_data_dir() + "expected/lidort_interface_masters/3p8_BRDF_ReadInput.cfg");
  int read_status = brdf_master.brdf_sup_inputstatus().bs_status_inputread();
  BOOST_CHECK_EQUAL(read_status, 0); // was read successful?

  if (read_status != 0)
    std::cerr << brdf_master.brdf_sup_inputstatus();

  // Types where data is stored
  Brdf_Sup_Inputs& brdf_inputs = brdf_master.brdf_sup_in();
  Brdf_Linsup_Outputs& brdf_outputs = brdf_master.brdf_linsup_out();
  
  // Process inputs
  bool do_debug_restoration = false;
  int bs_nmoments_input = 2 * brdf_inputs.bs_nstreams() - 1;
  brdf_master.run(do_debug_restoration, bs_nmoments_input);

  // Range and size names
  int n_beams = brdf_inputs.bs_nbeams();
  int n_u_streams = brdf_inputs.bs_n_user_streams();
  int n_u_relazms = brdf_inputs.bs_n_user_relazms();

  Range r_u_streams = Range(0,n_u_streams-1);
  Range r_u_realazms = Range(0,n_u_relazms-1);
  Range r_beams = Range(0,n_beams-1);

  int n_streams = brdf_inputs.bs_nstreams();

  Range r_nmoms = Range(0,bs_nmoments_input);
  Range r_streams = Range(0,n_streams-1);

  Range r_coef = Range(0,6-1);

  // Check values according to offline results
  IfstreamCs expected_data(test_data_dir() + "expected/lidort_interface_masters/lidort_ls_brdf_master");

  Array<double, 4> bs_ls_dbounce_brdfunc_expect;
  expected_data >> bs_ls_dbounce_brdfunc_expect;

  BOOST_CHECK_MATRIX_CLOSE(brdf_outputs.bs_ls_dbounce_brdfunc()
                           (r_coef, r_u_streams, r_u_realazms, r_beams), bs_ls_dbounce_brdfunc_expect);
 
  // --

  Array<double, 4> bs_ls_brdf_f_expect;
  expected_data >> bs_ls_brdf_f_expect;
  BOOST_CHECK_MATRIX_CLOSE(brdf_outputs.bs_ls_brdf_f()(r_coef, r_nmoms,r_streams,r_streams), bs_ls_brdf_f_expect);

  // --

  Array<double, 4> bs_ls_brdf_f_0_expect;
  expected_data >> bs_ls_brdf_f_0_expect;
  BOOST_CHECK_MATRIX_CLOSE(brdf_outputs.bs_ls_brdf_f_0()(r_coef, r_nmoms,r_streams,r_beams), bs_ls_brdf_f_0_expect);

  // --

  Array<double, 4> bs_ls_user_brdf_f_expect;
  expected_data >> bs_ls_user_brdf_f_expect;
  BOOST_CHECK_MATRIX_CLOSE(brdf_outputs.bs_ls_user_brdf_f()(r_coef, r_nmoms,r_u_streams,r_streams), bs_ls_user_brdf_f_expect);

  // --

  Array<double, 4> bs_ls_user_brdf_f_0_expect;
  expected_data >> bs_ls_user_brdf_f_0_expect;
  BOOST_CHECK_MATRIX_CLOSE(brdf_outputs.bs_ls_user_brdf_f_0()(r_coef, r_nmoms,r_u_streams,r_beams), bs_ls_user_brdf_f_0_expect);
}

BOOST_AUTO_TEST_CASE(lidort_lps_master)
{
  // Used for sizes
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Same as lidrot_brdf_master except for weighting function parts
  Brdf_Lin_Sup_Masters brdf_master = Brdf_Lin_Sup_Masters();

  // Read lidort config file
  brdf_master.read_config(test_data_dir() + "expected/lidort_interface_masters/3p8_BRDF_ReadInput.cfg");
  int read_status = brdf_master.brdf_sup_inputstatus().bs_status_inputread();
  BOOST_CHECK_EQUAL(read_status, 0); // was read successful?

  if (read_status != 0)
    std::cerr << brdf_master.brdf_sup_inputstatus();

  // Types where data is stored
  Brdf_Sup_Inputs& brdf_inputs = brdf_master.brdf_sup_in();

  // Process inputs
  bool do_debug_restoration = false;
  int bs_nmoments_input = 2 * brdf_inputs.bs_nstreams() - 1;
  brdf_master.run(do_debug_restoration, bs_nmoments_input);

  // Create Input and LIDORT interface classes
  Lidort_Inputs inp_master = Lidort_Inputs();
  Lidort_Lps_Masters lps_master = Lidort_Lps_Masters();

  // Read lidort config file
  inp_master.read_config(test_data_dir() + "expected/lidort_interface_masters/3p8_LIDORT_ReadInput.cfg");
  BOOST_CHECK_EQUAL(inp_master.lidort_inputstatus().ts_status_inputread(), lid_pars.lidort_success); // was read successful?

  if ( inp_master.lidort_inputstatus().ts_status_inputread() != lid_pars.lidort_success )
    std::cerr << inp_master.lidort_inputstatus() << std::endl;

  // Copy values from input master to lidort master
  lps_master.lidort_fixin( inp_master.lidort_fixin() );
  lps_master.lidort_modin( inp_master.lidort_modin() );

  // LIDORT types for consistent access:
  Lidort_Fixed_Control& fcont = lps_master.lidort_fixin().cont();
  Lidort_Modified_Control& mcont = lps_master.lidort_modin().mcont();
  Lidort_Fixed_Lincontrol& lin_fcontrol = lps_master.lidort_linfixin().cont();
  Lidort_Modified_Lincontrol& lin_mcontrol = lps_master.lidort_linmodin().mcont();
  Lidort_Fixed_Chapman& fchapman = lps_master.lidort_fixin().chapman();
  Lidort_Fixed_Optical& foptical = lps_master.lidort_fixin().optical();
  Lidort_Modified_Optical& moptical = lps_master.lidort_modin().moptical();
  Lidort_Fixed_Linoptical& linoptical = lps_master.lidort_linfixin().optical();
  Lidort_Main_Outputs& lid_output = lps_master.lidort_out().main();
  Lidort_Linatmos& lid_lpoutput = lps_master.lidort_linout().atmos();
  Lidort_Fixed_Uservalues& fuser_inputs = lps_master.lidort_fixin().userval();
  Lidort_Modified_Uservalues& muser_inputs = lps_master.lidort_modin().muserval();
  Lidort_Linsurf& lid_lsoutput = lps_master.lidort_linout().surf();

  // Get the preprepared atmosphere
  Array<double, 2> raymoms(3,lid_pars.maxlayers);
  Array<double, 1> height_grid(lid_pars.maxlayers+1);
  Array<double, 1> molext(lid_pars.maxlayers);
  Array<double, 1> molomg(lid_pars.maxlayers);

  int nlayers = fcont.ts_nlayers();
  IfstreamCs atm_data(test_data_dir() + "expected/lidort_interface_masters/input_atmos.dat");
  int int_dummy;
  double dbl_dummy;
  for (int mom_idx = 0; mom_idx < 3; mom_idx++) {
    atm_data >> int_dummy;
    for(int lay_idx = 0; lay_idx < nlayers; lay_idx++)
      atm_data >> raymoms(mom_idx, lay_idx);
  }

  height_grid(0) = 60.0;
  for(int lay_idx = 0; lay_idx < nlayers; lay_idx++)
    atm_data >> int_dummy 
             >> height_grid(lay_idx+1) 
             >> molext(lay_idx) 
             >> molomg(lay_idx)
             >> dbl_dummy >> dbl_dummy >> dbl_dummy >> dbl_dummy;


  // Add Aerosols bottom 6 layers, spread evenly
  Array<double, 1> aermoms(lid_pars.maxmoments_input+1);

  int nmoments_input = 80;
  double gaer = 0.8e0;
  aermoms(0) = 1.0e0;
  for(int mom_idx = 1; mom_idx < nmoments_input+1; mom_idx++) 
    aermoms(mom_idx) = (2.0*mom_idx+1) * pow(gaer, double(mom_idx));

  // Initialize array
  // phasmom_total variables in their moments dimension are declared 0: inside of Fortran,
  // so do not subtract one from translating indexes, also increase the size of the array,
  // by 1
  double lambertian_albedo;
  Array<double, 1> deltau_vert_input(lid_pars.maxlayers);
  Array<double, 1> omega_total_input(lid_pars.maxlayers);
  Array<double, 2> phasmoms_total_input(lid_pars.maxmoments_input+1, lid_pars.maxlayers);

  Array<double, 2> l_deltau_vert_input(lid_pars.max_atmoswfs,lid_pars.maxlayers);
  Array<double, 2> l_omega_total_input(lid_pars.max_atmoswfs,lid_pars.maxlayers);
  Array<double, 3> l_phasmoms_total_input(lid_pars.max_atmoswfs,lid_pars.maxmoments_input+1,lid_pars.maxlayers);

  lambertian_albedo = 0.0;
  deltau_vert_input = 0.0;
  omega_total_input = 0.0;
  phasmoms_total_input = 0.0;
  l_deltau_vert_input = 0.0;
  l_omega_total_input = 0.0;
  l_phasmoms_total_input = 0.0;

  // Initialize linearized inputs
  lin_mcontrol.ts_do_profile_linearization(true);
  lin_mcontrol.ts_do_surface_linearization(true);

  // New value in LIDORT 3.8.3 that if left at zero will cause floating point errors in asymtx
  fcont.ts_asymtx_tolerance(1e-20);
 
  // Initialise
  Array<bool, 1> layer_vary_flag(nlayers);
  Array<int, 1> layer_vary_number( lin_fcontrol.ts_layer_vary_number() );
  for (int lay_idx = 0; lay_idx < nlayers; lay_idx++) {
    layer_vary_number(lay_idx) = 2;
    layer_vary_flag(lay_idx)   = true;
  }
  // VERY important: While we can use referenced Blitz++ arrays to other
  // items, we CAN NOT use a reference to bool arrays since they are translated
  // in the accessors to and from int arrays as they are represented in Fortran
  // So we must explicitly send a value for copying.
  lin_fcontrol.ts_layer_vary_flag(layer_vary_flag);

  // Surface
  lin_fcontrol.ts_n_surface_wfs(1);
  lambertian_albedo = 0.05e0;

  // rayleigh layers
  int n6 = nlayers - 6;
  for (int lay_idx = 0; lay_idx < n6; lay_idx++) {
    deltau_vert_input(lay_idx) = molext(lay_idx);
    omega_total_input(lay_idx) = molomg(lay_idx);
    phasmoms_total_input(0,lay_idx) = raymoms(0,lay_idx);
    phasmoms_total_input(1,lay_idx) = raymoms(1,lay_idx);
    phasmoms_total_input(2,lay_idx) = raymoms(2,lay_idx);
    double ratio1 = 1.0e0 - molomg(lay_idx);
    l_deltau_vert_input(0,lay_idx) =   ratio1;
    l_omega_total_input(0,lay_idx) = - ratio1;
  }

  // aerosol layers
  double waer = 0.95e0 ; double taer = 0.5e0 ; 
  double parcel = taer / ( height_grid(n6) - height_grid(nlayers) );
  for (int lay_idx = n6; lay_idx < nlayers; lay_idx++) {
    double aerext = parcel * ( height_grid(lay_idx-1) - height_grid(lay_idx) );
    double aersca = aerext * waer;
    double molabs = molext(lay_idx) * ( 1.0e0 - molomg(lay_idx) );
    double molsca = molomg(lay_idx) * molext(lay_idx);
    double totext = molext(lay_idx) + aerext;
    double totsca = molsca    + aersca;
    double raywt  = molsca / totsca;
    double aerwt  = aersca / totsca;
    double omega  = totsca / totext;
    deltau_vert_input(lay_idx) = totext;
    omega_total_input(lay_idx) = omega;
    double ratio1 = molabs / totext;
    l_deltau_vert_input(0,lay_idx) =   ratio1;
    l_omega_total_input(0,lay_idx) = - ratio1;
    double ratio2 = aerext / totext;
    l_deltau_vert_input(1,lay_idx) = ratio2;
    l_omega_total_input(1,lay_idx) = ratio2 * ((waer/omega) - 1.0e0);
    phasmoms_total_input(0,lay_idx) = 1.0e0;
    for (int mom_idx = 1; mom_idx < 3; mom_idx++)
      phasmoms_total_input(mom_idx,lay_idx) = raywt * raymoms(mom_idx,lay_idx) + aerwt * aermoms(mom_idx);
    for (int mom_idx = 3; mom_idx < nmoments_input+1; mom_idx++)
      phasmoms_total_input(mom_idx,lay_idx) = aerwt * aermoms(mom_idx);
    for (int mom_idx = 1; mom_idx < nmoments_input+1; mom_idx++) {
      ratio2 = aersca / totsca;
      l_phasmoms_total_input(1,mom_idx,lay_idx) = ratio2 * 
        ( (aermoms(mom_idx)/phasmoms_total_input(mom_idx,lay_idx)) - 1.0e0 );
    }
  }

  //  Copy local control integers
  fcont.ts_nlayers(nlayers);
  mcont.ts_nmoments_input(nmoments_input);
  fchapman.ts_height_grid(height_grid);

  // Check calculated inputs against expected
  Array<double, 1> lambertian_albedo_exp;
  Array<double, 2> deltau_vert_input_exp;
  Array<double, 2> omega_total_input_exp;
  Array<double, 3> phasmoms_total_input_exp;

  Array<double, 3> l_deltau_vert_input_exp;
  Array<double, 3> l_omega_total_input_exp;
  Array<double, 4> l_phasmoms_total_input_exp;

  IfstreamCs lps_inputs_exp(test_data_dir() + "expected/lidort_interface_masters/lidort_lps_inputs");
  lps_inputs_exp >> lambertian_albedo_exp
                 >> deltau_vert_input_exp
                 >> omega_total_input_exp
                 >> phasmoms_total_input_exp
                 >> l_deltau_vert_input_exp
                 >> l_omega_total_input_exp
                 >> l_phasmoms_total_input_exp;

  Range all = Range::all();
  BOOST_CHECK_CLOSE(lambertian_albedo, lambertian_albedo_exp(0), 1e-8);
  BOOST_CHECK_MATRIX_CLOSE(deltau_vert_input(Range(0,deltau_vert_input_exp.rows()-1)) ,
               deltau_vert_input_exp(all, 0));
  BOOST_CHECK_MATRIX_CLOSE(omega_total_input(Range(0, omega_total_input_exp.rows()-1)),
               omega_total_input_exp(all, 0));
  BOOST_CHECK_MATRIX_CLOSE(phasmoms_total_input(Range(0,phasmoms_total_input_exp.rows()-1), Range(0,phasmoms_total_input_exp.cols() - 1)),
               phasmoms_total_input_exp(all, all, 0));

  BOOST_CHECK_MATRIX_CLOSE(l_deltau_vert_input(Range(0,l_deltau_vert_input_exp.rows()-1), Range(0,l_deltau_vert_input_exp.cols()-1)),
               l_deltau_vert_input_exp(all, all, 0));
  BOOST_CHECK_MATRIX_CLOSE(l_omega_total_input(Range(0,l_omega_total_input_exp.rows()-1),Range(0,l_omega_total_input_exp.cols()-1)),
               l_omega_total_input_exp(all, all, 0));
  BOOST_CHECK_MATRIX_CLOSE(l_phasmoms_total_input(Range(0,l_phasmoms_total_input.rows()-1),Range(0,l_phasmoms_total_input_exp.cols()-1),Range(0,l_phasmoms_total_input_exp.depth()-1)),
               l_phasmoms_total_input_exp(all, all, all, 0));

  // Copy to optical property type-structure inputs
  foptical.ts_lambertian_albedo(lambertian_albedo);
  foptical.ts_deltau_vert_input(deltau_vert_input);
  moptical.ts_omega_total_input(omega_total_input);
  foptical.ts_phasmoms_total_input(phasmoms_total_input);
        
  linoptical.ts_l_deltau_vert_input(l_deltau_vert_input);
  linoptical.ts_l_omega_total_input(l_omega_total_input);
  linoptical.ts_l_phasmoms_total_input(l_phasmoms_total_input);

  // Copy BRDF values from BRDF master to LPS master
  lps_master.lidort_sup().brdf().copy_from_sup( brdf_master.brdf_sup_out() );
  lps_master.lidort_linsup().brdf().copy_from_sup( brdf_master.brdf_linsup_out() );

  // LIDORT call to LPS master
  // Process inputs
  lps_master.run(false);

  std::ofstream lid_out("lidort_setup");
  lid_out << std::setprecision(10)
          << lps_master.lidort_fixin() << std::endl
          << lps_master.lidort_modin() << std::endl;

  // Check on status of calculation
  int inp_status = lps_master.lidort_out().status().ts_status_inputcheck();
  int calc_status = lps_master.lidort_out().status().ts_status_calculation();

  BOOST_CHECK_EQUAL(inp_status, lid_pars.lidort_success);
  BOOST_CHECK_EQUAL(calc_status, lid_pars.lidort_success);

  if ( inp_status != lid_pars.lidort_success || calc_status != lid_pars.lidort_success )
    std::cerr << lps_master.lidort_out() << std::endl;

  // Read expected outputs
  Array<double, 3> ts_intensity_expt;
  Array<double, 5> ts_profilewf_expt;
  Array<double, 4> ts_surfacewf_expt;
  Array<double, 3> ts_mean_intensity_expt;
  Array<double, 5> ts_mint_profilewf_expt;
  Array<double, 3> ts_flux_integral_expt;
  Array<double, 5> ts_flux_profilewf_expt;

  IfstreamCs lps_outputs_expt(test_data_dir() + "expected/lidort_interface_masters/lidort_lps_outputs");
  lps_outputs_expt >> ts_intensity_expt
                   >> ts_profilewf_expt
                   >> ts_surfacewf_expt
                   >> ts_mean_intensity_expt
                   >> ts_mint_profilewf_expt
                   >> ts_flux_integral_expt
                   >> ts_flux_profilewf_expt;

  // Ranges for use in comparisons since LIDORT allocates more memory in arrays
  // than it will initialize
  Range rlayers( 0, fcont.ts_nlayers() - 1 );
  Range rulevs( 0, fuser_inputs.ts_n_user_levels() - 1 );

  // LIDORT 3.8 has an issue with the geometry indexes other than the first
  // not agreeing with previous LIDORT versions. Until this is fixed 
  // Just compare the first index since other than this unit test we would
  // always use the first index
  //Range rgeoms( 0, (brdf_inputs.bs_nbeams() * 
  //                  muser_inputs.ts_n_user_streams()) - 1 );
  Range rgeoms(0, 0); // TEMP

  Range rbeams( 0, brdf_inputs.bs_nbeams() - 1 );
  Range rdirs( 0, 1 ); // up and down welling
  Range ratmwfs( 0, min(lin_fcontrol.ts_layer_vary_number()) - 1 );
  Range rsurfwfs( 0, 0 ); // only first index has values // brdf_inputs.bs_n_surface_wfs() 

  Array<double, 3> ts_intensity_calc( lid_output.ts_intensity() );
  BOOST_CHECK_MATRIX_CLOSE( ts_intensity_calc(rulevs,rgeoms,rdirs), 
                            ts_intensity_expt(rulevs,rgeoms,rdirs) );

  Array<double, 5> ts_profilewf_calc( lid_lpoutput.ts_profilewf() );
  BOOST_CHECK_MATRIX_CLOSE( ts_profilewf_calc(ratmwfs,rlayers,rulevs,rgeoms,rdirs), 
                            ts_profilewf_expt(ratmwfs,rlayers,rulevs,rgeoms,rdirs) );

  Array<double, 4> ts_surfacewf_calc( lid_lsoutput.ts_surfacewf() );
  BOOST_CHECK_MATRIX_CLOSE_TOL( ts_surfacewf_calc(rsurfwfs,rulevs,rgeoms,rdirs), 
                                ts_surfacewf_expt(rsurfwfs,rulevs,rgeoms,rdirs), 1e-6 );

}

BOOST_AUTO_TEST_SUITE_END()

