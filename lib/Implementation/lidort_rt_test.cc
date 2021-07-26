#include "lidort_rt.h"
#include "unit_test_support.h"
#include "lidort_fixture.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

void write_expected(std::string& filename, blitz::Array<double,1>& wn, ArrayAd<double,1>& refl_jac)
{
  std::ofstream outfile(filename.c_str());
  outfile << std::scientific
          << std::setprecision(20)
          << "# wavenumbers" << std::endl
          << wn << std::endl
          << "# reflectance" << std::endl
          << refl_jac.value() << std::endl
          << "# jacobian" << std::endl
          << refl_jac.jacobian() <<std::endl;
  outfile.close();
 
}

BOOST_FIXTURE_TEST_SUITE(lidort_rt_lambertian_compare, LidortLambertianFixture)

BOOST_AUTO_TEST_CASE(compare_sphericity)
{
  blitz::Array<double, 1> wn_expect;
  blitz::Array<double, 1> reflectance_expect;
  blitz::Array<double, 2> jacobian_expect;
  ArrayAd<double, 1> refl_jac_calc;

  // Turn off delta-m scaling
  lidort_rt->rt_driver()->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  ////////////////////// 
  // Plane parallel
  lidort_rt->rt_driver()->set_plane_parallel();

  refl_jac_calc.reference( lidort_rt->reflectance(wn_arr, 0).spectral_range().data_ad() );

  if(false) // Set to true to update expected results
    write_expected(pp_expected_filename, wn_arr, refl_jac_calc);
  IfstreamCs expected_data_pp(pp_expected_filename);
  expected_data_pp >> wn_expect >> reflectance_expect >> jacobian_expect;
  BOOST_CHECK_MATRIX_CLOSE(wn_expect, wn_arr);
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.value(), reflectance_expect, 1e-7);

  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.jacobian(), jacobian_expect, 1e-4);

  //////////////////////
  // Pseudo spherical
  lidort_rt->rt_driver()->set_pseudo_spherical();

  refl_jac_calc.reference( lidort_rt->reflectance(wn_arr, 0).spectral_range().data_ad() );

  if(false) // Set to true to update expected results
    write_expected(ps_expected_filename, wn_arr, refl_jac_calc);
  IfstreamCs expected_data_ps(ps_expected_filename);
  expected_data_ps >> wn_expect >> reflectance_expect >> jacobian_expect;

  BOOST_CHECK_MATRIX_CLOSE(wn_expect, wn_arr);
  BOOST_CHECK_MATRIX_CLOSE(refl_jac_calc.value(), reflectance_expect);

  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.jacobian(), jacobian_expect, 1e-4);
    
  //////////////////////
  // Plane parallel + sscorrection on
  lidort_rt->rt_driver()->set_plane_parallel_plus_ss_correction();

  refl_jac_calc.reference( lidort_rt->reflectance(wn_arr, 0).spectral_range().data_ad() );

  if(false) // Set to true to update expected results
    write_expected(pp_and_ss_expected_filename, wn_arr, refl_jac_calc);
  IfstreamCs expected_data_pp_and_ss(pp_and_ss_expected_filename);
  expected_data_pp_and_ss >> wn_expect >> reflectance_expect >> jacobian_expect;

  BOOST_CHECK_MATRIX_CLOSE(wn_expect, wn_arr);
  BOOST_CHECK_MATRIX_CLOSE(refl_jac_calc.value(), reflectance_expect);

  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.jacobian(), jacobian_expect, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE(lidort_rt_sphericity, LidortLambertianFixture)

BOOST_AUTO_TEST_CASE(check_los)
{
  int nstreams = 8;
  int nmoms = 180;
  bool do_multiple_scattering_only = false;

  Array<double, 1> vza;  
  Array<double, 1> stokes_ps;
  Array<double, 1> stokes_los;

  ArrayAd<double,2> empty_iv;

  // Fix sza, azm, but we will vary zen
  sza = 60;
  azm = 0.0;

  // From 0 to 89 degrees
  int num_angles = 24;

  vza.resize(num_angles);
  stokes_ps.resize(num_angles);
  stokes_los.resize(num_angles);
  for (int ang_idx = 0; ang_idx < num_angles; ang_idx++) {
    zen = 1e-6+ang_idx*(89.0/(num_angles-1));
    vza(ang_idx) = zen(0);

    if (zen(0) < 0.25) {
        pure_nadir = true;
    } else {
        pure_nadir = false;
    }

    // Set up PS driver
    boost::shared_ptr<LidortRt> lidort_ps;
    lidort_ps.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
                                 nstreams, nmoms, do_multiple_scattering_only));

    // Pseudo spherical
    lidort_ps->rt_driver()->set_pseudo_spherical();
    //lidort_ps->rt_driver()->lidort_interface()->lidort_modin().mbool().ts_do_no_azimuth(true);
    

    // Set up LOS driver, LOS should be default
    boost::shared_ptr<LidortRt> lidort_los;
    lidort_los.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
                                  nstreams, nmoms, do_multiple_scattering_only));

    // Run LIDORT drivers
    stokes_ps(ang_idx) = lidort_ps->stokes_single_wn(wn_arr(0), 0, empty_iv)(0);
    stokes_los(ang_idx) = lidort_los->stokes_single_wn(wn_arr(0), 0, empty_iv)(0);
  }

  // These two modes are similar up to a certain point
  Range similar_range = Range(0, last(vza < 60.0));
  BOOST_CHECK_MATRIX_CLOSE_TOL( stokes_ps(similar_range), stokes_los(similar_range), 1e-5 );

  // The point closest to 0 vza should be almost exactly the same
  BOOST_CHECK_CLOSE( stokes_ps(0), stokes_los(0), 1e-9 );

  /*
  std::ofstream outfile("los_compare.out");
  outfile << std::scientific
          << std::setprecision(20)
          << "# vza" << std::endl
          << vza << std::endl
          << "# ps" << std::endl
          << stokes_ps << std::endl
          << "# los" << std::endl
          << stokes_los <<std::endl;
  */
  // Python matplotlib code for plotting output above:
  // from matplotlib.pyplot import *
  // from Blitz_Array import Fast_Parser as Fpar
  // figure(1)
  // (vza, ps, los) = Fpar("los_compare.out").parse()
  // subplot(211);cla();plot(vza,ps);plot(vza,los);title("LIDORT 3.5 PS vs LOS - 8 streams, 180 moments");ylabel("Intensity");legend(["PS", "LOS"], "upper left")#;ylim([0.08,0.14])
  // subplot(212);cla();plot(vza,(ps-los)/los);title("LIDORT 3.5 (PS - LOS) / LOS for SZA = 60, AZM = 0");ylabel("Fractional difference");xlabel("Viewing Angle")#;ylim([0,0.0005])
  // show()
  //savefig("lidort_35_ps_los_compare.png")          

}

BOOST_AUTO_TEST_CASE(check_ps)
{
  int nstreams = 8;
  int nmoms = 180;
  bool do_multiple_scattering_only = false;

  Array<double, 1> chk_sza;  
  Array<double, 1> stokes_ps;
  Array<double, 1> stokes_pp;

  ArrayAd<double,2> empty_iv;

  // Fix zen, azm, but we will vary zen
  zen = 0.0001;
  azm = 0.0;

  // From 0 to 89 degrees
  int num_angles = 24;

  chk_sza.resize(num_angles);
  stokes_ps.resize(num_angles);
  stokes_pp.resize(num_angles);
  for (int ang_idx = 0; ang_idx < num_angles; ang_idx++) {
    sza = 1e-6+ang_idx*(89.0/(num_angles-1));
    chk_sza(ang_idx) = sza(0);

    // Set up PS driver
    boost::shared_ptr<LidortRt> lidort_ps;
    lidort_ps.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
                                 nstreams, nmoms, do_multiple_scattering_only));

    // Pseudo spherical
    lidort_ps->rt_driver()->set_pseudo_spherical();

    // Set up PP driver
    boost::shared_ptr<LidortRt> lidort_pp;
    lidort_pp.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
                                 nstreams, nmoms, do_multiple_scattering_only));

    // Plane-parallel + ss correction
    lidort_pp->rt_driver()->set_plane_parallel_plus_ss_correction();

    // Run LIDORT drivers
    stokes_ps(ang_idx) = lidort_ps->stokes_single_wn(wn_arr(0), 0, empty_iv)(0);
    stokes_pp(ang_idx) = lidort_pp->stokes_single_wn(wn_arr(0), 0, empty_iv)(0);
  }

  // These two modes are similar up to a certain point
  Range similar_range = Range(0, last(chk_sza < 40.0));
  BOOST_CHECK_MATRIX_CLOSE_TOL( stokes_ps(similar_range), stokes_pp(similar_range), 1e-4 );

  // The point closest to 0 sza should be almost exactly the same
  BOOST_CHECK_CLOSE( stokes_ps(0), stokes_pp(0), 1e-11 );

  /*
  std::ofstream outfile("ps_compare.out");
  outfile << std::scientific
          << std::setprecision(20)
          << "# chk_sza" << std::endl
          << chk_sza << std::endl
          << "# ps" << std::endl
          << stokes_ps << std::endl
          << "# pp" << std::endl
          << stokes_pp <<std::endl;
  */
  // Python matplotlib code for plotting output above:
  // from matplotlib.pyplot import *
  // from Blitz_Array import Fast_Parser as Fpar
  // figure(2)
  // (vza, ps, pp) = Fpar("ps_compare.out").parse()
  // subplot(211);cla();plot(vza,ps);plot(vza,pp);title("LIDORT 3.5 PS vs PP -  8 streams, 180 moments");ylabel("Intensity");legend(["PS", "PP"], "lower left");#ylim([0.08,0.14]);
  // subplot(212);cla();plot(vza,(ps-pp)/ps);title("LIDORT 3.5 (PS - PP) / PP for VZA = 0.0001, AZM = 0");ylabel("Fractional difference");xlabel("Solar Angle");#ylim([0,0.0005])
  // savefig("/home/mcduffie/Desktop/lidort_35_ps_pp_compare.png")
  // show()

}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////////////////////////////////////////////////////
// VVVVVVVVVVVVVV Cox-Munk VVVVVVVVVVVVVVV
///////////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE(lidort_rt_coxmunk_compare, LidortCoxmunkFixture)

BOOST_AUTO_TEST_CASE(compare_sphericity)
{
  blitz::Array<double, 1> wn_expect;
  blitz::Array<double, 1> reflectance_expect;
  blitz::Array<double, 2> jacobian_expect;

  ArrayAd<double, 1> refl_jac_calc;

  // Turn off delta-m scaling
  lidort_rt->rt_driver()->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  ////////////////////// 
  // Plane parallel
  lidort_rt->rt_driver()->set_plane_parallel();

  refl_jac_calc.reference( lidort_rt->reflectance(wn_arr, 0).spectral_range().data_ad() );

  IfstreamCs expected_data_pp(pp_expected_filename);
  expected_data_pp >> wn_expect >> reflectance_expect >> jacobian_expect;

  BOOST_CHECK_MATRIX_CLOSE(wn_expect, wn_arr);
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.value(), reflectance_expect,1e-6);

  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.jacobian(), jacobian_expect, 1e-7);

  //////////////////////
  // Pseudo spherical
  lidort_rt->rt_driver()->set_pseudo_spherical();

  refl_jac_calc.reference( lidort_rt->reflectance(wn_arr, 0).spectral_range().data_ad() );

  IfstreamCs expected_data_ps(ps_expected_filename);
  expected_data_ps >> wn_expect >> reflectance_expect >> jacobian_expect;

  BOOST_CHECK_MATRIX_CLOSE(wn_expect, wn_arr);
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.value(), reflectance_expect, 1e-6);

  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.jacobian(), jacobian_expect, 1e-6);

  //////////////////////
  // Plane parallel + sscorrection on
  lidort_rt->rt_driver()->set_plane_parallel_plus_ss_correction();

  refl_jac_calc.reference( lidort_rt->reflectance(wn_arr, 0).spectral_range().data_ad() );

  IfstreamCs expected_data_pp_and_ss(pp_and_ss_expected_filename);
  expected_data_pp_and_ss >> wn_expect >> reflectance_expect >> jacobian_expect;

  BOOST_CHECK_MATRIX_CLOSE(wn_expect, wn_arr);
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.value(), reflectance_expect, 1e-6);

  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.jacobian(), jacobian_expect, 1e-6);

}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////////////////////////////////////////////////////
// VVVVVVVVVVVVVV Cox-Munk Plus Lambertian VVVVVVVVVVVVVVV
///////////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE(lidort_rt_coxmunk_plus_lamb_compare, LidortCoxmunkPlusLambertianFixture)

BOOST_AUTO_TEST_CASE(compare_sphericity)
{
  blitz::Array<double, 1> wn_expect;
  blitz::Array<double, 1> reflectance_expect;
  blitz::Array<double, 2> jacobian_expect;

  ArrayAd<double, 1> refl_jac_calc;

  // Turn off delta-m scaling
  lidort_rt->rt_driver()->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);
  
  ////////////////////// 
  // Plane parallel
  lidort_rt->rt_driver()->set_plane_parallel();

  refl_jac_calc.reference( lidort_rt->reflectance(wn_arr, 0).spectral_range().data_ad() );

  IfstreamCs expected_data_pp(pp_expected_filename);
  expected_data_pp >> wn_expect >> reflectance_expect >> jacobian_expect;

  BOOST_CHECK_MATRIX_CLOSE(wn_expect, wn_arr);
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.value(), reflectance_expect, 1e-6);

  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.jacobian(), jacobian_expect, 1e-4);
  
  //////////////////////
  // Pseudo spherical
  lidort_rt->rt_driver()->set_pseudo_spherical();

  refl_jac_calc.reference( lidort_rt->reflectance(wn_arr, 0).spectral_range().data_ad() );

  IfstreamCs expected_data_ps(ps_expected_filename);
  expected_data_ps >> wn_expect >> reflectance_expect >> jacobian_expect;

  BOOST_CHECK_MATRIX_CLOSE(wn_expect, wn_arr);
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.value(), reflectance_expect, 1e-6);

  // Jacobians wrong for Coxmunk+Lambertian in LIDORT 3.0
  // 3.5 jacobians verified using Finite Difference test
  //BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.jacobian(), jacobian_expect, 1e-6);
      
  //////////////////////
  // Plane parallel + sscorrection on
  lidort_rt->rt_driver()->set_plane_parallel_plus_ss_correction();

  refl_jac_calc.reference( lidort_rt->reflectance(wn_arr, 0).spectral_range().data_ad() );

  IfstreamCs expected_data_pp_and_ss(pp_and_ss_expected_filename);
  expected_data_pp_and_ss >> wn_expect >> reflectance_expect >> jacobian_expect;

  BOOST_CHECK_MATRIX_CLOSE(wn_expect, wn_arr);
  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.value(), reflectance_expect, 1e-6);

  // Jacobians wrong for Coxmunk+Lambertian in LIDORT 3.0
  // 3.5 jacobians verified using Finite Difference test
  //BOOST_CHECK_MATRIX_CLOSE_TOL(refl_jac_calc.jacobian(), jacobian_expect, 1e-6);
  
}

BOOST_AUTO_TEST_SUITE_END()


