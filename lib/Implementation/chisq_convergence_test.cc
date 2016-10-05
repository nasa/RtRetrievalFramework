#include "chisq_convergence.h"
#include "unit_test_support.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(chisq_convergence, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  ChisqConvergence conv;
  FitStatistic fs1, fs2;
  bool has_converged, convergence_failed, step_diverged;
  double gamma = 0.0;
  fs2.number_iteration = 1;
  conv.convergence_check(fs1, fs2, has_converged, convergence_failed, 
			 gamma, step_diverged);
  BOOST_CHECK_EQUAL(has_converged, false);
  BOOST_CHECK_EQUAL(convergence_failed, false);
  BOOST_CHECK_EQUAL(step_diverged, false);
  BOOST_CHECK_CLOSE(1.0 + gamma, 1.0, 1e-8);
  fs2.chisq_measured = 0.1;
  fs1.chisq_measured = 1.0;
  fs2.number_iteration = 2;
  conv.convergence_check(fs1, fs2, has_converged, convergence_failed, 
			 gamma, step_diverged);
  BOOST_CHECK_EQUAL(has_converged, false);
  BOOST_CHECK_EQUAL(convergence_failed, false);
  BOOST_CHECK_EQUAL(step_diverged, false);
  BOOST_CHECK_CLOSE(1.0 + gamma, 0.1, 1e-8);
  gamma = 0.0;
  fs2.chisq_measured = 1.0;
  fs1.chisq_measured = 0.1;
  conv.convergence_check(fs1, fs2, has_converged, convergence_failed, 
			 gamma, step_diverged);
  BOOST_CHECK_EQUAL(has_converged, false);
  BOOST_CHECK_EQUAL(convergence_failed, false);
  BOOST_CHECK_EQUAL(step_diverged, true);
  BOOST_CHECK_CLOSE(1.0 + gamma, 10.0, 1e-8);
  gamma = 0.0;
  fs2.number_iteration = 50;
  conv.convergence_check(fs1, fs2, has_converged, convergence_failed, 
			 gamma, step_diverged);
  BOOST_CHECK_EQUAL(has_converged, false);
  BOOST_CHECK_EQUAL(convergence_failed, true);
  BOOST_CHECK_EQUAL(step_diverged, false);
  BOOST_CHECK_CLOSE(1.0 + gamma, 1.0, 1e-8);
  fs2.chisq_measured = 1e-8;
  fs2.number_iteration = 2;
  conv.convergence_check(fs1, fs2, has_converged, convergence_failed, 
			 gamma, step_diverged);
  BOOST_CHECK_EQUAL(has_converged, true);
  BOOST_CHECK_EQUAL(convergence_failed, false);
  BOOST_CHECK_EQUAL(step_diverged, false);
  BOOST_CHECK_CLOSE(1.0 + gamma, 1.0, 1e-8);
  fs2.chisq_measured = fs1.chisq_measured * 0.99999999;
  conv.convergence_check(fs1, fs2, has_converged, convergence_failed, 
			 gamma, step_diverged);
  BOOST_CHECK_EQUAL(has_converged, true);
  BOOST_CHECK_EQUAL(convergence_failed, false);
  BOOST_CHECK_EQUAL(step_diverged, false);
  BOOST_CHECK_CLOSE(1.0 + gamma, 1.0, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
