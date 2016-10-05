#include "composite_initial_guess.h"
#include "unit_test_support.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(composite_initial_guess, GlobalFixture)

class InitialGuessBuilderTest : public InitialGuessBuilder {
public:
  InitialGuessBuilderTest(int start, int num) :s(start), n(num) {}
  virtual ~InitialGuessBuilderTest() {}
  virtual int number_element() const {return n;}
  virtual void build_initial_value(blitz::Array<double,1>& v, int index) const
  {
    for(int i = 0; i < n; ++i)
      v(index + i) = s + i;
  }
  virtual void build_apriori(blitz::Array<double,1>& v, int index) const
  {
    for(int i = 0; i < n; ++i)
      v(index + i) = s + 2 * i;
  }
  virtual void build_apriori_covariance(blitz::Array<double, 2>& m, 
					int index) const
  {
    for(int i = 0; i < n; ++i)
      m(index + i, index + i) = s + i;
  }
private:
  int s, n;
};
BOOST_AUTO_TEST_CASE(basic)
{
  CompositeInitialGuess ci;
  boost::shared_ptr<InitialGuessBuilder> i1(new InitialGuessBuilderTest(0,3));
  boost::shared_ptr<InitialGuessBuilder> i2(new InitialGuessBuilderTest(7,4));
  ci.add_builder(i1);
  ci.add_builder(i2);
  blitz::Array<double, 1> expect(7);
  blitz::Array<double, 2> expect2(7,7);
  expect = 0,1,2,7,8,9,10;
  BOOST_CHECK_MATRIX_CLOSE(ci.initial_guess(), expect);
  expect = 0,2,4,7,9,11,13;
  BOOST_CHECK_MATRIX_CLOSE(ci.apriori(), expect);
  expect2 =
    0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,
    0,0,2,0,0,0,0,
    0,0,0,7,0,0,0,
    0,0,0,0,8,0,0,
    0,0,0,0,0,9,0,
    0,0,0,0,0,0,10;
  BOOST_CHECK_MATRIX_CLOSE(ci.apriori_covariance(), expect2);
  ci.remove_builder(i1);
  expect.resize(4);
  expect2.resize(4,4);
  expect = 7,8,9,10;
  BOOST_CHECK_MATRIX_CLOSE(ci.initial_guess(), expect);
  expect = 7,9,11,13;
  BOOST_CHECK_MATRIX_CLOSE(ci.apriori(), expect);
  expect2 =
    7,0,0,0,
    0,8,0,0,
    0,0,9,0,
    0,0,0,10;
  BOOST_CHECK_MATRIX_CLOSE(ci.apriori_covariance(), expect2);
}

BOOST_AUTO_TEST_SUITE_END()
