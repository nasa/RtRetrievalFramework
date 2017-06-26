#include "composite_perturbation.h"
#include "unit_test_support.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(composite_perturbation, GlobalFixture)

class PerturbationBuilderTest : public PerturbationBuilder {
public:
  PerturbationBuilderTest(int start, int num) :s(start), n(num) {}
  virtual ~PerturbationBuilderTest() {}
  virtual int number_element() const {return n;}
  virtual void build_perturbation(blitz::Array<double,1>& v, int index) const
  {
    for(int i = 0; i < n; ++i)
      v(index + i) = s + i;
  }
private:
  int s, n;
};
BOOST_AUTO_TEST_CASE(basic)
{
  CompositePerturbation ci;
  boost::shared_ptr<PerturbationBuilder> i1(new PerturbationBuilderTest(0,3));
  boost::shared_ptr<PerturbationBuilder> i2(new PerturbationBuilderTest(7,4));
  ci.add_builder(i1);
  ci.add_builder(i2);
  blitz::Array<double, 1> expect(7);
  expect = 0,1,2,7,8,9,10;
  BOOST_CHECK_MATRIX_CLOSE(ci.perturbation(), expect);
  ci.remove_builder(i1);
  expect.resize(4);
  expect = 7,8,9,10;
  BOOST_CHECK_MATRIX_CLOSE(ci.perturbation(), expect);
}

BOOST_AUTO_TEST_SUITE_END()
