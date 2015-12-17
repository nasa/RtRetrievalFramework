#include "l2_fp_configuration_lua.h"
#include "radiative_transfer.h"
#include "output_hdf.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(rayleigh_only, GlobalFixture)

BOOST_AUTO_TEST_CASE(radiance_and_jacobian)
{
  // Compare a Rayleigh only radiative transfer with one that has
  // aerosols in it, but all the extinction coefficients set to zero.
  L2FpConfigurationLua c_zeroext(test_data_dir() + "config_zeroext.lua");
  L2FpConfigurationLua c_rayleigh(test_data_dir() + "config_rayleigh_only.lua");
  const RadiativeTransfer& rt_zeroext = 
    *c_zeroext.lua_state().globals()["config"]["rt"].
    value_ptr<RadiativeTransfer>();
  const RadiativeTransfer& rt_rayleigh = 
    *c_rayleigh.lua_state().globals()["config"]["rt"].
    value_ptr<RadiativeTransfer>();
  Array<double, 1> wn_arr(10);
  wn_arr = 12930, 12940, 12950, 12960, 12970, 12980, 12990, 
    13000, 13010, 13020;
  ArrayAd<double, 1> refl_zeroext = 
    rt_zeroext.reflectance(wn_arr, 0).spectral_range().data_ad();
  ArrayAd<double, 1> refl_rayleigh = 
    rt_rayleigh.reflectance(wn_arr, 0).spectral_range().data_ad();
  BOOST_CHECK_MATRIX_CLOSE(refl_zeroext.value(), refl_rayleigh.value());
  Array<double, 2> jac_zeroext = refl_zeroext.jacobian();
  Array<double, 2> jac_rayleigh = refl_rayleigh.jacobian();
  // We need to subset jac_zeroext to take out the aerosol stuff in
  // order to match up to jac_rayleigh.
  Array<double, 2> jac_zeroext_sub(jac_rayleigh.shape());
  Range ra(Range::all());
  jac_zeroext_sub(ra, Range(0, 22)) = jac_zeroext(ra, Range(0,22));
  jac_zeroext_sub(ra, Range(23,31)) = jac_zeroext(ra, Range(111-(31-23),111));
  BOOST_CHECK_MATRIX_CLOSE(jac_zeroext_sub, jac_rayleigh);

  // Check that we also can write out the results.
  boost::shared_ptr<Output> out(new OutputHdf("rayleigh_only_output.h5", 20, 
					      jac_rayleigh.cols(), 5, 3));
  boost::shared_ptr<Output> outerr(new OutputHdf("rayleigh_only_output.h5.error", 20, 
						 jac_rayleigh.cols(), 5, 3));
  add_file_to_cleanup("rayleigh_only_output.h5");
  c_rayleigh.output(out, outerr);
  IfstreamCs in(test_data_dir() + "connor_rayleigh_only.txt");
  in >> *(c_rayleigh.solver());
  c_rayleigh.forward_model()->state_vector()->update_state(c_rayleigh.solver()->x_solution(), 
					  c_rayleigh.solver()->aposteriori_covariance());
  out->write();
}

BOOST_AUTO_TEST_SUITE_END()
