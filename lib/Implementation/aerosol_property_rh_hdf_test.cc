#include "aerosol_property_rh_hdf.h"
#include "atmosphere_fixture.h"
#include "pressure_sigma.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_property_rh_hdf, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile h(test_data_dir() + "l2_merra_aerosol_RH.h5");
  AerosolPropertyRhHdf a(h, "SS/Properties", 
			 atm->pressure_ptr(), atm->relative_humidity_ptr());
  BOOST_CHECK_CLOSE(a.extinction_coefficient_each_layer(13000).value()(0), 
		    93.839991056848092, 1e-8);
  BOOST_CHECK_CLOSE(a.extinction_coefficient_each_layer(13000).value()(17), 
		    144.7355511600864, 1e-8);
  BOOST_CHECK_CLOSE(a.scattering_coefficient_each_layer(13000).value()(0), 
		    93.828162933185922, 1e-8);
  BOOST_CHECK_CLOSE(a.scattering_coefficient_each_layer(13000).value()(17), 
		    144.7286068121517, 1e-8);
  // Phase function is large, so we just check the first couple of moments.
  Array<double, 2> pf_expect(2, 6);
  pf_expect =
    1, 0.00536656, 0.911908, 0.0145849, 1, 0.911908,
    2.37414, 0.0742056, 2.4446, -0.0642392, 2.37414, 2.4446;
  if(false)
    std::cerr << a.phase_function_moment_each_layer(13000).value()
      (Range(0,1), 0, Range::all()) << "\n";
  BOOST_CHECK_MATRIX_CLOSE_TOL(a.phase_function_moment_each_layer(13000).value()
			       (Range(0,1), 0, Range::all()), pf_expect, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()

