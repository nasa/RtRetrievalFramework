#include "initial_guess_value.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(initial_guess_value, GlobalFixture)
BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<InitialGuessValue> ig(new InitialGuessValue());
  Array<double, 1> v(2);
  v = 1, 2;
  Array<double, 2> m(2,2);
  m = 1, 2, 3, 4;
  ig->apriori(v);
  ig->apriori_covariance(m);
  BOOST_CHECK_MATRIX_CLOSE(ig->apriori(), v);
  BOOST_CHECK_MATRIX_CLOSE(ig->initial_guess(), v);
  BOOST_CHECK_MATRIX_CLOSE(ig->apriori_covariance(), m);
  boost::shared_ptr<InitialGuessValue> ig2(new InitialGuessValue());
  Array<double, 1> v2(3);
  v2 = 1, 2, 3;
  Array<double, 2> m2(3,3);
  m2 = 1, 2, 3, 
    4, 5, 6,
    7, 8, 9;
  ig2->apriori(v2);
  ig2->apriori_covariance(m2);
  CompositeInitialGuess ci;
  ci.add_builder(ig);
  ci.add_builder(ig2);
  blitz::Array<double, 1> expect(5);
  blitz::Array<double, 2> expect2(5,5);
  expect = 1,2,1,2,3;
  BOOST_CHECK_EQUAL(ci.number_element(), 5);
  BOOST_CHECK_MATRIX_CLOSE(ci.initial_guess(), expect);
  BOOST_CHECK_MATRIX_CLOSE(ci.apriori(), expect);
  expect2 = 
    1, 2, 0, 0, 0,
    3, 4, 0, 0, 0,
    0, 0, 1, 2, 3,
    0, 0, 4, 5, 6,
    0, 0, 7, 8, 9;
  BOOST_CHECK_MATRIX_CLOSE(ci.apriori_covariance(), expect2);
  v = 4, 5;
  ig->initial_guess(v);
  BOOST_CHECK_MATRIX_CLOSE(ci.apriori(), expect);
  expect = 4, 5, 1, 2, 3;
  BOOST_CHECK_MATRIX_CLOSE(ci.initial_guess(), expect);
}

BOOST_AUTO_TEST_SUITE_END()
