#include "composite_initial_guess.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(composite_initial_guess_full_guess, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(full_initial_guess)
{
  const InitialGuess& ig = *config_initial_guess;
  IfstreamCs in(test_data_dir() + 
		"expected/composite_initial_guess/full_initial_guess");
  Array<double, 1> initial_sv_expected, apriori_sv_expected;
  Array<double, 2> apriori_cov_expected;
  in >> initial_sv_expected >> apriori_sv_expected >> apriori_cov_expected;
  BOOST_CHECK_MATRIX_CLOSE(ig.initial_guess(), initial_sv_expected);
  BOOST_CHECK_MATRIX_CLOSE(ig.apriori(), apriori_sv_expected);
  BOOST_CHECK_MATRIX_CLOSE(ig.apriori_covariance(), apriori_cov_expected);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(composite_initial_guess_full_guess_ecmwf, ConfigurationEcmwfFixture)

BOOST_AUTO_TEST_CASE(full_initial_guess_ecmwf)
{
  const InitialGuess& ig = *config_initial_guess;
  IfstreamCs in(test_data_dir() + 
		"expected/composite_initial_guess/full_initial_guess_ecmwf");
  Array<double, 1> initial_sv_expected, apriori_sv_expected;
  Array<double, 2> apriori_cov_expected;
  in >> initial_sv_expected >> apriori_sv_expected >> apriori_cov_expected;
  BOOST_CHECK_MATRIX_CLOSE(ig.initial_guess(), initial_sv_expected);
  BOOST_CHECK_MATRIX_CLOSE(ig.apriori(), apriori_sv_expected);
  BOOST_CHECK_MATRIX_CLOSE(ig.apriori_covariance(), apriori_cov_expected);
}

BOOST_AUTO_TEST_SUITE_END()

