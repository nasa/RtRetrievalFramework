#include "aerosol_shape_gaussian.h"
#include "aerosol_shape_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_shape_gaussian, AerosolShapeFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<Pressure> p = config_pressure;

  // Arrays defining coeffs
  blitz::Array<bool, 1> ret_flag(3);
  ret_flag = true;
  blitz::Array<double, 1> coeffs(3);

  // Loaded expected results
  IfstreamCs gaussian_expt(test_data_dir() + "expected/aerosol_shape/gaussian");

  // Kahn
  Array<double, 1> kahn_expt;
  gaussian_expt >> kahn_expt;

  coeffs = -3.28341, 1.0, 0.2;
  AerosolShapeGaussian aer_kahn_log = AerosolShapeGaussian(p, ret_flag, coeffs, "Kahn", false);
  BOOST_CHECK_MATRIX_CLOSE_TOL(kahn_expt, aer_kahn_log.aerosol_extinction().value(), 1e-14);

  coeffs(0) = exp(coeffs(0));
  AerosolShapeGaussian aer_kahn_lin = AerosolShapeGaussian(p, ret_flag, coeffs, "Kahn", true);
  BOOST_CHECK_MATRIX_CLOSE_TOL(kahn_expt, aer_kahn_lin.aerosol_extinction().value(), 1e-14);

  // Water
  Array<double, 1> water_expt;
  gaussian_expt >> water_expt;

  coeffs = -3.28341, 0.75, 0.1;
  AerosolShapeGaussian aer_water_log = AerosolShapeGaussian(p, ret_flag, coeffs, "Water", false);
  BOOST_CHECK_MATRIX_CLOSE_TOL(water_expt, aer_water_log.aerosol_extinction().value(), 1e-14);

  coeffs(0) = exp(coeffs(0));
  AerosolShapeGaussian aer_water_lin = AerosolShapeGaussian(p, ret_flag, coeffs, "Water", true);
  BOOST_CHECK_MATRIX_CLOSE_TOL(water_expt, aer_water_lin.aerosol_extinction().value(), 1e-14);

  // Ice
  Array<double, 1> ice_expt;
  gaussian_expt >> ice_expt;

  coeffs = -3.28341, 0.3, 0.04;
  AerosolShapeGaussian aer_ice_log = AerosolShapeGaussian(p, ret_flag, coeffs, "Ice", false);
  BOOST_CHECK_MATRIX_CLOSE_TOL(ice_expt, aer_ice_log.aerosol_extinction().value(), 1e-14);

  coeffs(0) = exp(coeffs(0));
  AerosolShapeGaussian aer_ice_lin = AerosolShapeGaussian(p, ret_flag, coeffs, "Ice", true);
  BOOST_CHECK_MATRIX_CLOSE_TOL(ice_expt, aer_ice_lin.aerosol_extinction().value(), 1e-14);
}

BOOST_AUTO_TEST_CASE(config)
{
  boost::shared_ptr<AerosolOptical> aer = 
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);

  // Loaded expected results
  IfstreamCs gaussian_expt(test_data_dir() + "expected/aerosol_shape/gaussian");
  IfstreamCs gaussian_jac_expt(test_data_dir() + "expected/aerosol_shape/gaussian_jacobian");

  std::vector<std::string> aer_name = aer->aerosol_name();
  for(int aer_idx = 1; aer_idx < aer->number_particle(); aer_idx++) {
    boost::shared_ptr<AerosolExtinctionImpBase> aer_obj(boost::dynamic_pointer_cast<AerosolExtinctionImpBase>(aer->aerosol_extinction(aer_idx)));
    ArrayAd<double, 1> aext = aer_obj->aerosol_extinction();

    Array<double, 1> value_expt;
    gaussian_expt >> value_expt;
    Array<double, 2> jac_expt;
    gaussian_jac_expt >> jac_expt;
    BOOST_CHECK_MATRIX_CLOSE_TOL(value_expt, aext.value(), 1e-14);
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac_expt, aext.jacobian(), 1e-14);  
  }
}

BOOST_AUTO_TEST_SUITE_END()

