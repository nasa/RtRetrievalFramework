#include "twostream_rt.h"
#include "lidort_fixture.h"
#include "unit_test_support.h"
#include "fe_disable_exception.h"

using namespace FullPhysics;
using namespace blitz;

void compare_lidort_2stream(boost::shared_ptr<RtAtmosphere>& atmosphere, 
                  const boost::shared_ptr<StokesCoefficient>& stokes_coefs,
			    blitz::Array<double, 1>& sza, 
                            blitz::Array<double, 1>& zen, 
                            blitz::Array<double, 1>& azm, 
                            blitz::Array<double, 1>& wn_arr,
                            bool debug_output)
{ 
  int nstreams = 1;
  int nmoms = 3;
  bool do_multiple_scattering_only = true;
  bool pure_nadir = false;

  boost::shared_ptr<LidortRt> lidort_rt;
  lidort_rt.reset(new LidortRt(atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
                               nstreams, nmoms, do_multiple_scattering_only));
  
  // Turn off delta-m scaling
  lidort_rt->rt_driver()->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  lidort_rt->rt_driver()->set_plane_parallel();

  ArrayAd<double, 2> lid_rad_jac = lidort_rt->stokes_and_jacobian(wn_arr, 0);

  boost::shared_ptr<TwostreamRt> twostream_rt;
  twostream_rt.reset(new TwostreamRt(atmosphere, stokes_coefs, sza, zen, azm, false));

  // Turn off delta-m scaling
  twostream_rt->rt_driver()->twostream_interface()->do_d2s_scaling(false);

  // Plane-parallel
  twostream_rt->rt_driver()->twostream_interface()->do_plane_parallel(true);

  ArrayAd<double, 2> ts_rad_jac = twostream_rt->stokes_and_jacobian(wn_arr, 0);

  Range all = Range::all();
  if(debug_output) {
    std::cerr << "lid_rad = " << lid_rad_jac.value()(all, 0) << std::endl
              << "ts_rad = " << ts_rad_jac.value()(all, 0) << std::endl;
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(lid_rad_jac.value()(all, 0), ts_rad_jac.value()(all, 0), 1e-8);

  for(int jac_idx = 0; jac_idx < ts_rad_jac.jacobian().depth(); jac_idx++) {
    if(debug_output) {
      std::cerr << jac_idx << " -- l: ";
      for(int wn_idx = 0; wn_idx < lid_rad_jac.jacobian().rows(); wn_idx++)      
        std::cerr << lid_rad_jac.jacobian()(wn_idx, 0, jac_idx) << " ";
      std::cerr << std::endl << jac_idx << " -- t: ";
      for(int wn_idx = 0; wn_idx < ts_rad_jac.jacobian().rows(); wn_idx++)      
        std::cerr << ts_rad_jac.jacobian()(wn_idx, 0, jac_idx) << " ";
      std::cerr << std::endl;
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(lid_rad_jac.jacobian()(all, 0, jac_idx), ts_rad_jac.jacobian()(all, 0, jac_idx), 3e-5);
  }
}

BOOST_FIXTURE_TEST_SUITE(twostream_rt_lambertian, LidortLambertianFixture)

BOOST_AUTO_TEST_CASE(comparison)
{ 
  compare_lidort_2stream(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, false);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(twostream_rt_coxmunk, LidortCoxmunkFixture)

BOOST_AUTO_TEST_CASE(comparison)
{ 
  compare_lidort_2stream(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, false);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(twostream_rt_coxmunk_plus_lambertian, LidortCoxmunkPlusLambertianFixture)

BOOST_AUTO_TEST_CASE(comparison)
{ 
  // Not sure what is triggering an exception on ubuntu. It is hard to 
  // duplicate this problem, so we'll just shut this off for now.
  FeDisableException disable_fp;
  compare_lidort_2stream(config_atmosphere, stokes_coefs, sza, zen, azm, wn_arr, false);
}

BOOST_AUTO_TEST_SUITE_END()

