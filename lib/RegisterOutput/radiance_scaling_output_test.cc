#include "unit_test_support.h"
#include "radiance_scaling_output.h"
#include "output_hdf.h"
#include "configuration_fixture.h"
#include "radiance_scaling_sv_fit.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(radiance_scaling_output, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  blitz::Array<double, 1> coeff(2);
  coeff = 0.5, 0.0;
  blitz::Array<bool, 1> flag(2);
  flag = true;
  DoubleWithUnit band_ref(1e4/0.77, "cm^-1");
  std::string band_name = "o2";
  boost::shared_ptr<RadianceScaling> rad_scaling(new RadianceScalingSvFit(coeff, flag, band_ref, band_name));

  RadianceScalingOutput po(rad_scaling, band_name);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("radiance_scaling_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("radiance_scaling_output.h5");
  po.register_output_apriori(out);
  po.register_output(out);

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()
