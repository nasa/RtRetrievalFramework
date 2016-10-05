#include "lsi_rt.h"
#include "fp_logger.h"
#include "lidort_fixture.h"
#include "unit_test_support.h"
#include "hdf_file.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(lsi_rt, LidortLowHighLambertianFixture)

BOOST_AUTO_TEST_CASE(err_est)
{
  is_long_test();		// Skip unless we are running long tests.
  turn_on_logger();		// Have log output show up.
  
  HdfFile config(test_data_dir() + "l2_fixed_level_static_input.h5");
  IfstreamCs expected_data(test_data_dir() + "expected/lsi_rt/lsi_expected");
  LsiRt rt(low_rt, high_rt, config);
  int wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01)
    wn_i += 1;
  Array<double, 1> wn_arr(wn_i);
  wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01, ++wn_i)
    wn_arr(wn_i) = wn;
  ArrayAd<double, 2> err_est = rt.correction_only(wn_arr, 0);
  Logger::info() << low_rt->atmosphere_ptr()->timer_info();
  Array<double, 1> err_est_expect;
  expected_data >> err_est_expect;
  if(false) {			// Write to output, if we need to
				// regenerate expected values.
    std::cerr.precision(20);
    std::cerr << err_est.value()(Range::all(), 0) << "\n";
  }
  BOOST_CHECK_MATRIX_CLOSE_TOL(err_est.value()(Range::all(), 0), 
			       err_est_expect, 1e-7);
}

BOOST_AUTO_TEST_CASE(err_est_jac)
{
  // Even for a long test, this takes a long time to run (several
  // minutes). We don't normally run this, since all the jacobian
  // calculation are done by AutoDerivative anyways which
  // automatically calculates this. But we can turn this test on if
  // needed.
  is_really_long_test();     // Skip unless we are running long tests.

  HdfFile c(test_data_dir() + "l2_fixed_level_static_input.h5");
  IfstreamCs expected_data(test_data_dir() + "expected/lsi_rt/lsi_expected");
  LsiRt rt(low_rt, high_rt, c);

  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(config_initial_guess->initial_guess());
  sv.update_state(sv0);

  int spec_index = 2;
  int wn_i = 0;
  for(double wn = 4810; wn <= 4897; wn += 0.01)
    wn_i += 1;
  Array<double, 1> wn_arr(wn_i);
  wn_i = 0;
  for(double wn = 4810; wn <= 4897; wn += 0.01, ++wn_i)
    wn_arr(wn_i) = wn;

  ArrayAd<double, 2> err_est = rt.correction_only(wn_arr, spec_index);
  Array<double, 2> e0(err_est.shape());
  e0 = err_est.value();
  Array<double, 3> jac = err_est.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(err_est.shape());
    jacfd = (rt.correction_only(wn_arr, spec_index).value() - e0) 
      / epsilon(i);
    Array<double, 2> diff(jac(Range::all(), Range::all(), i) - jacfd);
    if(false) {			// Can turn this on to dump values,
				// if needed for debugging
      if(max(abs(diff)) !=0) {
	std::cerr << i << ": " << max(abs(diff)) << " "
		  << max(abs(where(abs(jacfd) < 1e-15, 0, diff/jacfd))) << "\n";
      }
    }
    if(i < 20)
    // The CO2 VMR jacobian values are much larger than the rest, so
    // tolerance if more reasonable for those sizes
      BOOST_CHECK(max(abs(diff)) < 0.2);
    else
      BOOST_CHECK(max(abs(diff)) < 3e-5);
  }
}

BOOST_AUTO_TEST_CASE(stokes)
{
  // We don't normally need to run this test. All the real calculation
  // is in err_est test above. Stokes is just a direct function of the
  // err_est. We need a test to make sure that the actual calculation
  // is ok, but given the time it takes to run this we don't need to
  // generally do so.
  is_really_long_test();     // Skip unless we are running long tests.
  turn_on_logger();		// Have log output show up.
  
  HdfFile config(test_data_dir() + "l2_fixed_level_static_input.h5");
  IfstreamCs expected_data(test_data_dir() + "expected/lsi_rt/stokes");
  LsiRt rt(low_rt, high_rt, config);
  int wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01)
    wn_i += 1;
  Array<double, 1> wn_arr(wn_i);
  wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01, ++wn_i)
    wn_arr(wn_i) = wn;
  Array<double, 1> refl = rt.reflectance(wn_arr, 0, true).spectral_range().data();
  Array<double, 1> refl_expect;
  expected_data >> refl_expect;
  BOOST_CHECK_MATRIX_CLOSE(refl(Range(9,19)), refl_expect);
  expected_data >> refl_expect;
  BOOST_CHECK_MATRIX_CLOSE(refl(Range(15999,16009)), refl_expect);
  if(false) {			// Write to output, if we need to
				// regenerate expected values.
    std::cerr.precision(20);
    std::cerr << "# refl Range(9,19)" << std::endl
	      << std::scientific << refl(Range(9,19)) << std::endl
	      << "# refl Range(15999,16009)" << std::endl
	      << std::scientific << refl(Range(15999,16009)) << std::endl;
  }

}

BOOST_AUTO_TEST_CASE(stokes_and_jacobian)
{
  // We don't normally need to run this test. All the real calculation
  // is in err_est test above. Stokes is just a direct function of the
  // err_est. We need a test to make sure that the actual calculation
  // is ok, but given the time it takes to run this we don't need to
  // generally do so.
  is_really_long_test();     // Skip unless we are running long tests.
  turn_on_logger();		// Have log output show up.

  IfstreamCs expected_data(test_data_dir() + 
			   "expected/lsi_rt/stokes");
  std::string lsi_fname = test_data_dir() + "old_ascii/lsi_wl.dat";
  LsiRt rt(low_rt, high_rt, lsi_fname);
  int wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01)
    wn_i += 1;
  Array<double, 1> wn_arr(wn_i);
  wn_i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01, ++wn_i)
    wn_arr(wn_i) = wn;
  ArrayAd<double, 1> refl = rt.reflectance(wn_arr, 0).spectral_range().data_ad();
  Array<double, 1> refl_expect;
  expected_data >> refl_expect;
  BOOST_CHECK_MATRIX_CLOSE(refl.value()(Range(9,19)), refl_expect);
  expected_data >> refl_expect;
  BOOST_CHECK_MATRIX_CLOSE(refl.value()(Range(15999,16009)), refl_expect);
  if(false) {			// Write to output, if we need to
				// regenerate expected values.
    std::cerr.precision(20);
    std::cerr << "# refl Range(9,19)" << std::endl
	      << std::scientific << refl.value()(Range(9,19)) << std::endl
	      << "# refl Range(15999,16009)" << std::endl
	      << std::scientific << refl.value()(Range(15999,16009)) << std::endl;
  }

}

BOOST_AUTO_TEST_SUITE_END()

// ConfigurationFixture that reads the input file 
// config_20100207003330.lua. This is a test case extracted from real
// data that lets us test the LSI for the all 3 stokes parameters (our
// standard case above has coefficient 1,0,0 so it ignores the 2nd and
// 3rd stokes coefficients).

class Configuration20100207003330 : public ConfigurationFixture {
public:
  Configuration20100207003330() 
    : ConfigurationFixture("config_20100207003330.lua") {}
  virtual ~Configuration20100207003330() {}
};

BOOST_FIXTURE_TEST_SUITE(lsi_rt2, Configuration20100207003330)
BOOST_AUTO_TEST_CASE(lsi_all_stokes)
{
  is_long_test();              // Skip unless we are running long tests.
  turn_on_logger();            // Have log output show up.
  int band = 1;		       // Look at weak CO2 (this had a problem
			       // described in ticket #946)
  SpectralDomain spec_domain = highres_grid(band);
  ArrayAd<double, 1> refl = config_rt->reflectance(spec_domain,band).spectral_range().data_ad();
  // Look at temperature offset. This was value initial plotted in
  // problem described in ticket #946. This still is a good one to
  // look at.
  int jac_index = 22;
  if(false) {			// If we need to regenerate output.
    std::cerr.precision(12);
    std::cerr << refl.value() << "\n";
    std::cerr << Array<double, 1>(refl.jacobian()(Range::all(), jac_index));
  }
  IfstreamCs expected(test_data_dir() + "expected/lsi_rt/lsi_all_stokes");
  Array<double, 1> refl_expect;
  Array<double, 1> jac_expect;
  expected >> refl_expect >> jac_expect;
  BOOST_CHECK_MATRIX_CLOSE(refl.value(), refl_expect);
  BOOST_CHECK_MATRIX_CLOSE(refl.jacobian()(Range::all(), jac_index),
			   jac_expect);
}

// This next text is in response to Ticket #939. When we run the LSI
// twice on the same input, we expect to get the same results. For the
// very first point (and only for the very first point), we had a
// condition where the Jacobians were difference, primarily for the
// CO2 VMR. This test illustrates this problem.
BOOST_AUTO_TEST_CASE(lsi_run_twice)
{
  is_long_test();              // Skip unless we are running long tests.
  turn_on_logger();            // Have log output show up.
  int band = 1;		       // Look at weak CO2 
  SpectralDomain spec_domain = highres_grid(band);
  ArrayAd<double, 1> refl = config_rt->reflectance(spec_domain,band).spectral_range().data_ad();
  ArrayAd<double, 1> refl2 = config_rt->reflectance(spec_domain,band).spectral_range().data_ad();
  BOOST_CHECK_MATRIX_CLOSE(refl.jacobian(), refl2.jacobian());
}
BOOST_AUTO_TEST_SUITE_END()
