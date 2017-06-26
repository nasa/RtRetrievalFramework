#include "temperature_met.h"
#include "atmosphere_fixture.h"
#include "acos_sounding_id.h"
#include "unit_test_support.h"
#include "temperature_level_offset.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(temperature_level_offset, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::string sid = "20091009203401";
  HdfFile sfile(test_data_dir() + "in/sounding_id.h5");
  std::vector<boost::shared_ptr<HdfSoundingId> > sidv = 
    AcosSoundingId::create(sfile, sid);
  boost::shared_ptr<Pressure> p = config_pressure;
  StateVector sv;

  Array<double, 1> temp_expect(19);
  temp_expect = 244.2, 214.553, 218.029, 222.544, 218.341, 221.37, 227.38,
    233.493, 239.376, 244.52, 248.708, 251.979, 254.537, 256.655, 258.521,
    260.155, 261.747, 261.732, 258.598;

  TemperatureLevelOffset t1(p, temp_expect, 0, true);
  sv.add_observer(t1);

  for(int i = 0; i < temp_expect.rows(); ++i) {
    BOOST_CHECK_CLOSE(t1.temperature(p->pressure_grid()(i)).convert(units::K).value.value(), temp_expect(i), 1e-3);
  }
}

BOOST_AUTO_TEST_SUITE_END()

