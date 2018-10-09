#include "unit_test_support.h"
#include "temperature_met_output.h"
#include "output_hdf.h"
#include "acos_sounding_id.h"
#include "configuration_fixture.h"
#include "acos_met_file.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(temperature_met_output, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::string sid = "20091009203401";
  HdfFile sfile(test_data_dir() + "in/sounding_id.h5");
  std::vector<boost::shared_ptr<HdfSoundingId> > sidv = 
    AcosSoundingId::create(sfile, sid);
  boost::shared_ptr<AcosMetFile> e(new AcosMetFile(test_data_dir() + "in/ecmwf.h5", 
					       sidv[0], sidv.size() > 1));
  boost::shared_ptr<Pressure> p = config_pressure;
  TemperatureMetOutput po(boost::shared_ptr<TemperatureMet>(new TemperatureMet(e, p, 0, true)));
  boost::shared_ptr<OutputHdf> out(new OutputHdf("temperature_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("temperature_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in temperature unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


