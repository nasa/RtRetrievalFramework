#include "aerosol_property_hdf.h"
#include "pressure_sigma.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_property_hdf, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile h(test_data_dir() + "l2_fixed_level_static_input.h5");
  Array<double, 1> a1(3), b(3);
  a1 = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  boost::shared_ptr<Pressure> p(new PressureSigma(a1,b, psurf, true));
  AerosolPropertyHdf a(h, "Aerosol/Kahn_2b/Properties", p);
  BOOST_CHECK_CLOSE(a.extinction_coefficient_each_layer(13000).value()(0), 
		    0.9321898305, 1e-8);
  BOOST_CHECK_CLOSE(a.scattering_coefficient_each_layer(13000).value()(0), 
		    0.877556390083, 1e-8);
  // Phase function is large, so we just check the first couple of moments.
  Array<double, 2> pf_expect(2, 6);
  pf_expect =
    1, 0, 0, 0.895197915, 0, 0,
    2.178746875, 0, 0, 2.24237365833, 0, 0;
  BOOST_CHECK_MATRIX_CLOSE(a.phase_function_moment_each_layer(13000).value()
			   (Range(0,1), 0, Range::all()), pf_expect);
}

BOOST_AUTO_TEST_SUITE_END()

