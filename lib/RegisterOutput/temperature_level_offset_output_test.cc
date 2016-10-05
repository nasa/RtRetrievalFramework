#include "unit_test_support.h"
#include "temperature_level_offset_output.h"
#include "output_hdf.h"
#include "acos_sounding_id.h"
#include "configuration_fixture.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(temperature_level_offset_output, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::string sid = "20091009203401";
  HdfFile sfile(test_data_dir() + "in/sounding_id.h5");
  std::vector<boost::shared_ptr<HdfSoundingId> > sidv = 
    AcosSoundingId::create(sfile, sid);
  boost::shared_ptr<Pressure> p = config_pressure;

  Array<double, 1> t1(19);
  t1 = 244.2, 214.553, 218.029, 222.544, 218.341, 221.37, 227.38,
    233.493, 239.376, 244.52, 248.708, 251.979, 254.537, 256.655, 258.521,
    260.155, 261.747, 261.732, 258.598;
  
  TemperatureLevelOffsetOutput po(boost::shared_ptr<TemperatureLevelOffset>(new TemperatureLevelOffset(p, t1, 0, true)));
  boost::shared_ptr<OutputHdf> out(new OutputHdf("temperature_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("temperature_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in temperature unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


