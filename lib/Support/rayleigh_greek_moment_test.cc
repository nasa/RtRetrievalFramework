#include "rayleigh_greek_moment.h"
#include "unit_test_support.h"
#include "hdf_file.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(rayleigh_greek_moment, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile h(test_data_dir() + "l2_fixed_level_static_input.h5");
  double depolar_fact = h.read_attribute<double>("/Gas/Air/depolarization_factor");
  blitz::Array<double, 2> array_expect(3,6);
  array_expect = 
    1,     0, 0, 0, 0, 0,
    1e-11, 0, 0, 1.3968144385817844, 0, 0,
    0.47936288771635682, 2.8761773262981407, 0, 0, 1.1741944765321404, 0;
  BOOST_CHECK(max(abs(RayleighGreekMoment::array(depolar_fact) - array_expect)) < 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
