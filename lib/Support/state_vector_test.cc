#include "state_vector.h"
#include "fp_exception.h"
#include "unit_test_support.h"
#include <iostream>
using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(state_vector, GlobalFixture)

class StateVectorTest : public SubStateVectorObserver {
public:
  StateVectorTest()
    : SubStateVectorObserver(3)
  {
  }
  virtual void mark_used_sub(blitz::Array<bool, 1>& Used) const
  {
    Used(1) = true; Used(2) = true;
  }
  virtual ~StateVectorTest() {}
  virtual void update_sub_state
  (const ArrayAd<double, 1>& Sv_sub,
   const blitz::Array<double, 2>& Cov)
  {
    x.reference(Sv_sub.value());
  }
  blitz::Array<double, 1> x;
  virtual void print(std::ostream& Os) const { Os << "StateVectorTest"; }
};

BOOST_AUTO_TEST_CASE(basic)
{
  StateVector sv;
  StateVectorTest t1;
  sv.add_observer(t1);
  StateVectorTest t2;
  sv.add_observer(t2);
  blitz::Array<double, 1> x(6);
  x = 1, 2, 3, 4, 5, 6;
  sv.update_state(x);
  BOOST_CHECK_EQUAL(sv.state_vector_name()(0), "State vector 1");
  BOOST_CHECK_EQUAL(sv.state_vector_name()(1), "State vector 2");
  BOOST_CHECK_EQUAL(sv.state_vector_name()(2), "State vector 3");
  BOOST_CHECK_EQUAL(sv.state_vector_name()(3), "State vector 4");
  BOOST_CHECK_EQUAL(sv.state_vector_name()(4), "State vector 5");
  BOOST_CHECK_EQUAL(sv.state_vector_name()(5), "State vector 6");
  blitz::Array<double, 1> xexpect(3);
  xexpect = 1, 2, 3;
  BOOST_CHECK_MATRIX_CLOSE(t1.x, xexpect);
  xexpect = 4, 5, 6;
  BOOST_CHECK_MATRIX_CLOSE(t2.x, xexpect);
  BOOST_CHECK_MATRIX_CLOSE(sv.state(), x);
  BOOST_CHECK_EQUAL(sv.used_flag()(0), false);
  BOOST_CHECK_EQUAL(sv.used_flag()(1), true);
  BOOST_CHECK_EQUAL(sv.used_flag()(2), true);
  BOOST_CHECK_EQUAL(sv.used_flag()(3), false);
  BOOST_CHECK_EQUAL(sv.used_flag()(4), true);
  BOOST_CHECK_EQUAL(sv.used_flag()(5), true);
}

BOOST_AUTO_TEST_CASE(bad_data)
{
  StateVector sv;
  StateVectorTest t1;
  sv.add_observer(t1);
  StateVectorTest t2;
  sv.add_observer(t2);
  blitz::Array<double, 1> xbad(5);
  BOOST_CHECK_THROW(sv.update_state(xbad), Exception);
  blitz::Array<double, 1> x(7);
  x = 1, 2, 3, 4, 5, 6, 7;
  sv.update_state(x);
}

BOOST_AUTO_TEST_SUITE_END()
