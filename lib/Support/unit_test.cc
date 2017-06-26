#include "unit.h"
#include "fp_exception.h"
#include "unit_test_support.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(unit, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  BOOST_CHECK_EQUAL(units::m.name(), "m");
  BOOST_CHECK_CLOSE(FullPhysics::conversion(units::m, units::cm), 100.0, 1e-8);
  BOOST_CHECK_THROW(FullPhysics::conversion(units::kg, units::m), Exception);
}

BOOST_AUTO_TEST_CASE(parser)
{
  using namespace units;
  // Handle prefix
  BOOST_CHECK(Unit::parse("ym").is_commensurate(m));
  BOOST_CHECK_CLOSE(Unit::parse("ym").conversion_to_si(), 1e-24,
		    1e-8);
  // Handle base type without a prefix
  BOOST_CHECK(Unit::parse("m").is_commensurate(m));
  BOOST_CHECK_CLOSE(Unit::parse("m").conversion_to_si(), 1.0,
		    1e-8);
  // Handle multiple with spaces
  BOOST_CHECK(Unit::parse("m s s").is_commensurate(m * s * s));
  // Handle () and "*"
  BOOST_CHECK(Unit::parse("m (s * s)").is_commensurate(m * s * s));
  // Handle ^
  BOOST_CHECK(Unit::parse("m s ^ 2").is_commensurate(m * s * s));
  // And negative power
  BOOST_CHECK(Unit::parse("m s ^ -2").is_commensurate(m / (s * s)));
  // And grouping
  BOOST_CHECK(Unit::parse("m (kg s)^2").
	      is_commensurate(m * kg * kg * s * s));
  // Don't allow base units together without a " " or "*"
  BOOST_CHECK_THROW(Unit::parse("ss"), Exception);
  // Handle "/"
  BOOST_CHECK(Unit::parse("m / (kg s)^2").
	      is_commensurate(m / kg / kg / s / s));
  // Handle combination of /.
  BOOST_CHECK(Unit::parse("m / kg / s^2").
	      is_commensurate(m / kg / (s * s)));
  // Handle combination of / with grouping
  BOOST_CHECK(Unit::parse("m / (kg / s^2)").
	      is_commensurate(m / kg * (s * s)));
}

BOOST_AUTO_TEST_CASE(equality_test)
{
  Unit u1("m / (kg s)^2");
  Unit u2("m / kg^2 / s^2");
  Unit u3("m / g^2 / s^2");
  Unit u4("m / kg / s^2");
  BOOST_CHECK_EQUAL(u1, u2);
  BOOST_CHECK(u1 != u3);
  BOOST_CHECK(u1 != u4);
}

BOOST_AUTO_TEST_SUITE_END()

