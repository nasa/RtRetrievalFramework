#include "unit_test_support.h"
#include "oco_sounding_id.h"
#include "oco_sim_met_ecmwf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(oco_sim_met_ecmwf, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::string sid_s = "2010090912004075";
  HdfFile h(test_data_dir() + "oco2_sim_l1b.h5");
  boost::shared_ptr<OcoSoundingId> sid(new OcoSoundingId(h, sid_s));
  OcoSimMetEcmwf f(test_data_dir() + "oco2_sim_met.h5", sid);
  Array<double, 1> press_grid = f.pressure_levels();
  Array<double, 1> temp_grid = f.temperature();
  BOOST_CHECK_CLOSE(press_grid(10), 177.81173706054688, 1e-4);
  BOOST_CHECK_CLOSE(temp_grid(10), 216.83796691894531,  1e-4);
  Array<double, 1> spec_hum = f.specific_humidity();
  BOOST_CHECK_CLOSE(spec_hum(10), 2.6958341550198384e-06,  1e-4);
  BOOST_CHECK_CLOSE(f.surface_pressure(), 86876.6640625, 1e-4);
  BOOST_CHECK_CLOSE(f.windspeed(), 2.1595962047576904, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()

