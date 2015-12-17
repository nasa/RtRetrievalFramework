#include "aerosol_property_hdf.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_property_hdf, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile h(test_data_dir() + "l2_fixed_level_static_input.h5");
  AerosolPropertyHdf a(h, "Aerosol/Kahn_2b/Properties");
  BOOST_CHECK_CLOSE(a.extinction_coefficient(13000), 0.9321898305, 1e-8);
  BOOST_CHECK_CLOSE(a.scattering_coefficient(13000), 0.877556390083, 1e-8);
  // Phase function is large, so we just check the first couple of moments.
  Array<double, 2> pf_expect(2, 6);
  pf_expect =
    1, 0, 0, 0.895197915, 0, 0,
    2.178746875, 0, 0, 2.24237365833, 0, 0;
  BOOST_CHECK_MATRIX_CLOSE(a.phase_function_moment(13000)(Range(0,1), 
							  Range::all()),
			   pf_expect);
}

BOOST_AUTO_TEST_SUITE_END()

