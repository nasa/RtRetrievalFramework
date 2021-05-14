#include "unit_test_support.h"
#include "l_rad_rt.h"
#include "lsi_rt.h"
#include "lidort_fixture.h"

using namespace FullPhysics;
using namespace blitz;

// Perhaps move this elsewhere
class LRadLambertianFixture : public LidortLambertianFixture {
public:
  LRadLambertianFixture()
  {
    // Set viewing geometry
    sza = 74.128288268999995;
    zen = 30.0;
    azm = 10.0;

    rt_first_order.reset(new LRadRt(lidort_rt, config_spectral_window->spectral_bound(), 
                                    sza, zen, azm, pure_nadir, true, false));
    rt_second_order.reset(new LRadRt(lidort_rt, config_spectral_window->spectral_bound(),
                                     sza, zen, azm, pure_nadir, true, true));
    rt_lrad_only.reset(new LRadRt(lidort_rt->stokes_coefficient(),
                                  config_atmosphere,
                                  config_spectral_window->spectral_bound(),
                                  sza, zen, azm, pure_nadir,
                                  lidort_rt->number_stokes()));
    rt_lrad_only_second.reset(new LRadRt(lidort_rt->stokes_coefficient(),
                                  config_atmosphere,
                                  config_spectral_window->spectral_bound(),
                                  sza, zen, azm, pure_nadir,
                                  lidort_rt->number_stokes(), true));
    check_tol = 1e-6;
  }
  virtual ~LRadLambertianFixture() {}

  double check_tol;
  boost::shared_ptr<LRadRt> rt_first_order;
  boost::shared_ptr<LRadRt> rt_second_order;
  boost::shared_ptr<LRadRt> rt_lrad_only;
  boost::shared_ptr<LRadRt> rt_lrad_only_second;
};

class LRadCoxmunkFixture : public LidortCoxmunkFixture {
public:
  LRadCoxmunkFixture()
  {
    // Set viewing geometry
    sza = 74.128288268999995;
    zen = 30.0;
    azm = 10.0;

    rt_first_order.reset(new LRadRt(lidort_rt, config_spectral_window->spectral_bound(), 
                                        sza, zen, azm, pure_nadir, true, false));
    rt_second_order.reset(new LRadRt(lidort_rt, config_spectral_window->spectral_bound(),
                                         sza, zen, azm, pure_nadir, true, true));
    rt_lrad_only.reset(new LRadRt(lidort_rt->stokes_coefficient(),
                                      config_atmosphere,
                                      config_spectral_window->spectral_bound(),
                                      sza, zen, azm, pure_nadir,
                                      lidort_rt->number_stokes()));
    rt_lrad_only_second.reset(new LRadRt(lidort_rt->stokes_coefficient(),
                                      config_atmosphere,
                                      config_spectral_window->spectral_bound(),
                                      sza, zen, azm, pure_nadir,
                                      lidort_rt->number_stokes(), true));
    check_tol = 1e-6;
  }
  virtual ~LRadCoxmunkFixture() {}

  double check_tol;
  boost::shared_ptr<LRadRt> rt_first_order;
  boost::shared_ptr<LRadRt> rt_second_order;
  boost::shared_ptr<LRadRt> rt_lrad_only;
  boost::shared_ptr<LRadRt> rt_lrad_only_second;
};

// --------- lambertian ----------- //

BOOST_FIXTURE_TEST_SUITE(l_rad_rt_lambertian, LRadLambertianFixture)

BOOST_AUTO_TEST_CASE(reflectance_first_order)
{
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.077138484532828524576, -0.00095639739803483338979, -0.0003442980540008311937, 
    0.077136521983769981703, -0.00095650473289548632646, -0.00034433669398848179826;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order->stokes(wn, 0),
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order->stokes_and_jacobian(wn, 0).value(),
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_lrad_only)
{
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.028119871962910821667, -0.00095639739803483338979, -0.0003442980540008311937,
    0.028119063326642042971, -0.00095650473289548632646, -0.00034433669398848179826;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_lrad_only->stokes(wn, 0),
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_lrad_only->stokes_and_jacobian(wn, 0).value(),
                               stokes_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_second_order)
{
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.077130330079821715628, -0.001102530809783065259, -0.00040364944623159254057,
    0.077128365808671203729, -0.0011026544975300596897, -0.00040369468472197715756;

  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order->stokes(wn, 0), 
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order->stokes_and_jacobian(wn, 0).value(), 
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(l_rad_timing)
{
  is_timing_test();
  boost::timer tm;
  RtAtmosphere& atm = *config_atmosphere;
  int i = 0;
  ArrayAd<double, 2> dummy(0,0);
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01) {
    ArrayAd<double, 1> 
      res(rt_lrad_only->stokes_and_jacobian_single_wn(wn, 0, dummy));
    if(++i % 1000 == 0)
      std::cerr << "Done with " << i << "\n"
                << "Total: " << tm.elapsed() << "\n"
                << atm.timer_info() << "\n";
  }
}

BOOST_AUTO_TEST_CASE(jac_first_order)
{
  // We use rt_lrad_only so we only test the portion of the Jacobian coming
  // from LRad. We separately test LIDORT.
  is_long_test();
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(config_initial_guess->initial_guess());
  sv.update_state(sv0);
  int spec_index = 2;
  Array<double, 1> wn(1);
  wn = 4820.0;

  // Test jacobians for all three l_rad_first pseudo spherical modes
  for(int ps_idx = LRadDriver::REGULAR; ps_idx <= LRadDriver::PLANE_PARALLEL; ps_idx++) {
    LRadDriver::PsMode ps_mode = static_cast<LRadDriver::PsMode>(ps_idx);
    boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                         config_atmosphere,
                                                         config_spectral_window->spectral_bound(),
                                                         sza, zen, azm, pure_nadir,
                                                         rt_lrad_only->number_stokes(),
                                                         false, 4, 0.01,
                                                         ps_mode));
    ArrayAd<double, 2> stk = 
      lrad_ps->stokes_and_jacobian(wn, spec_index);
    Array<double, 2> stk0(stk.shape());
    stk0 = stk.value();
    Array<double, 3> jac = stk.jacobian().copy();
    for(int i = 0; i < sv.state().rows(); ++i) {
      Array<double, 1> svn(sv0.copy());
      svn(i) += epsilon(i);
      sv.update_state(svn);
      Array<double, 2> jacfd(stk.shape());
      jacfd = (lrad_ps->stokes(wn, spec_index) - stk0) / epsilon(i);
      double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
      if(false) {                        // Can turn this on to dump values,
                                  // if needed for debugging
        if(diff > 0) {
          std::cerr << i << ": " << diff << "\n"
                    << diff / max(abs(jacfd)) << "\n"
                    << jac(0,Range::all(), i) << "\n"
                    << jacfd << "\n";
        }
      }
      if(diff > 1e-6) {
        std::cerr << i << ": " << diff << "\n"
                  << diff / max(abs(jacfd)) << "\n"
                  << jac(0,Range::all(), i) << "\n"
                  << jacfd << "\n";
      }
      BOOST_CHECK(diff < 1e-6);
    }
  }
}

BOOST_AUTO_TEST_CASE(jac_second_order)
{
  // We use rt_lrad_only so we only test the portion of the Jacobian coming
  // from LRad. We separately test LIDORT.
  is_long_test();
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(config_initial_guess->initial_guess());
  sv.update_state(sv0);
  int spec_index = 2;
  Array<double, 1> wn(1);
  wn = 4820.0;
  ArrayAd<double, 2> stk = 
    rt_lrad_only_second->stokes_and_jacobian(wn, spec_index);
  Array<double, 2> stk0(stk.shape());
  stk0 = stk.value();
  Array<double, 3> jac = stk.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(stk.shape());
    jacfd = (rt_lrad_only_second->stokes(wn, spec_index) - stk0) 
      / epsilon(i);
    double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
    if(false) {                        // Can turn this on to dump values,
                                // if needed for debugging
      if(diff > 0) {
        std::cerr << i << ": " << diff << "\n"
                  << diff / max(abs(jacfd)) << "\n"
                  << jac(0,Range::all(), i) << "\n"
                  << jacfd << "\n";
      }
    }
    if(diff > 1e-4) {
      std::cerr << i << ": " << diff << "\n"
                << diff / max(abs(jacfd)) << "\n"
                << jac(0,Range::all(), i) << "\n"
                << jacfd << "\n";
    }
    BOOST_CHECK(diff < 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(regular_ps)
{
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir,
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::REGULAR));
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.027538673425655725568, -0.00088752226179599554446, -0.00032303168554235871116;
  
  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(plane_parallel_ps)
{
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir,
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::PLANE_PARALLEL));
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.027411149233159647076, -0.00088605720558977132603, -0.00032249844869176630575;

  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

// --------- coxmunk ----------- //

BOOST_FIXTURE_TEST_SUITE(l_rad_rt_coxmunk, LRadCoxmunkFixture)

BOOST_AUTO_TEST_CASE(reflectance_first_order)
{
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.0079613385061960122852, -0.0015461389239466109130, -0.00055660191445895516611,
    0.0079615098335930490486, -0.0015462247310868153412, -0.00055663280457999967689;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order->stokes(wn, 0), 
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order->stokes_and_jacobian(wn, 0).value(),
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_first_order->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_first_order->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(reflectance_second_order)
{
  Array<double, 1> wn(2);
  wn = 12929.94, 12930.30;
  Array<double, 2> stokes_expect(2, 3);
  stokes_expect = 
    0.0079681282901736590757, -0.0024250036984472403759, -0.00074540308734331189606,
    0.0079682994187089527943, -0.0024250909594306377728, -0.00074543862544189246414;
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order->stokes(wn, 0), 
                               stokes_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order->stokes_and_jacobian(wn, 0).value(), 
     stokes_expect, check_tol);
  Array<double, 1> rad_expect(stokes_expect(Range::all(), 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(rt_second_order->reflectance(wn, 0, true).
                               spectral_range().data(), 
                               rad_expect, check_tol);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (rt_second_order->reflectance(wn, 0).spectral_range().data(), 
     rad_expect, check_tol);
}

BOOST_AUTO_TEST_CASE(jac_first_order)
{
  // We use rt_lrad_only so we only test the portion of the Jacobian coming
  // from LRad. We separately test LIDORT.
  is_long_test();
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(config_initial_guess->initial_guess());
  sv.update_state(sv0);

  int spec_index = 2;
  Array<double, 1> wn(1);
  wn = 4820.0;

  // Test jacobians for all three l_rad_first pseudo spherical modes
  for(int ps_idx = LRadDriver::REGULAR; ps_idx <= LRadDriver::PLANE_PARALLEL; ps_idx++) {
    LRadDriver::PsMode ps_mode = static_cast<LRadDriver::PsMode>(ps_idx);
    boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                         config_atmosphere,
                                                         config_spectral_window->spectral_bound(),
                                                         sza, zen, azm, pure_nadir,
                                                         rt_lrad_only->number_stokes(),
                                                         false, 4, 0.01,
                                                         ps_mode));

    ArrayAd<double, 2> stk = 
      lrad_ps->stokes_and_jacobian(wn, spec_index);
    Array<double, 2> stk0(stk.shape());
    stk0 = stk.value();
    Array<double, 3> jac = stk.jacobian().copy();
    for(int i = 103; i < sv.state().rows(); ++i) {
    //  for(int i = 0; i < sv.state().rows(); ++i) {
      Array<double, 1> svn(sv0.copy());
      svn(i) += epsilon(i);
      sv.update_state(svn);
      Array<double, 2> jacfd(stk.shape());
      jacfd = (lrad_ps->stokes(wn, spec_index) - stk0) 
        / epsilon(i);
      double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
      // Can turn this on to dump values,
      // if needed for debugging
      bool debug = false;
      if(debug || (diff > 1e-6)) {
        std::cerr << "SV:" << i << ", Initial guess: " << sv0(i) << std::endl
                  << "epsilon = " << epsilon(i) << std::endl
                  << "abs diff = " << diff << std::endl;
        double fdmax_diff = max(abs(jacfd));
        if (fdmax_diff > 0.0)
          std::cerr << "rel diff = "<< diff / fdmax_diff << std::endl;
        std::cerr << "analytic: " << jac(0,Range::all(), i) << std::endl
                  << "finite diff: " << jacfd << std::endl;
      }
      BOOST_CHECK(diff < 1e-6);
    }
  }
}

BOOST_AUTO_TEST_CASE(jac_second_order)
{
  // We use rt_lrad_only so we only test the portion of the Jacobian coming
  // from LRad. We separately test LIDORT.
  is_long_test();
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(config_initial_guess->initial_guess());
  sv.update_state(sv0);
  int spec_index = 2;
  Array<double, 1> wn(1);
  wn = 4820.0;
  ArrayAd<double, 2> stk = 
    rt_lrad_only_second->stokes_and_jacobian(wn, spec_index);
  Array<double, 2> stk0(stk.shape());
  stk0 = stk.value();
  Array<double, 3> jac = stk.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(stk.shape());
    jacfd = (rt_lrad_only_second->stokes(wn, spec_index) - stk0) 
      / epsilon(i);
    double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
    // Can turn this on to dump values,
    // if needed for debugging
    bool debug = false;
    if(debug || (diff > 1e-4)) {
      std::cerr << "SV:" << i << ", Initial guess: " << sv0(i) << std::endl
                << "epsilon = " << epsilon(i) << std::endl
                << "abs diff = " << diff << std::endl;
      double fdmax_diff = max(abs(jacfd));
      if (fdmax_diff > 0.0)
        std::cerr << "rel diff = "<< diff / fdmax_diff << std::endl;
      std::cerr << "analytic: " << jac(0,Range::all(), i) << std::endl
                << "finite diff: " << jacfd << std::endl;
    }
    BOOST_CHECK(diff < 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(regular_ps)
{
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir,
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::REGULAR));
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.0023771646969834359965, -0.00088753741227118679042, -0.00032303719986433816411;

  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(plane_parallel_ps)
{
  // Set up to run l_rad in regular ps mode
  Array<double, 1> zen_small(3);
  zen_small = 1.0e-6;
  boost::shared_ptr<LRadRt> lrad_ps(new LRadRt(rt_lrad_only->stokes_coefficient(),
                                                        config_atmosphere,
                                                        config_spectral_window->spectral_bound(),
                                                        sza, zen_small, azm, pure_nadir, 
                                                        rt_lrad_only->number_stokes(),
                                                        false, 4, 0.01,
                                                        LRadDriver::PLANE_PARALLEL));
  Array<double, 1> wn(1);
  wn = 12929.94;
  Array<double, 1> stokes_expect(3);
  stokes_expect = 
    0.0023728677739882617949, -0.00088607228186708588812, -0.00032250393600792706195;

  BOOST_CHECK_MATRIX_CLOSE_TOL(lrad_ps->stokes(wn, 0)(0, Range::all()), stokes_expect, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
