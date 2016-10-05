#include "oco_sim_apriori.h"
#include "hdf_file.h"
#include "oco_sounding_id.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(oco_sim_apriori, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::string sid_s = "2010090912004075";
  HdfFile h(test_data_dir() + "oco2_sim_l1b.h5");
  OcoSoundingId sid(h, sid_s);
  OcoSimApriori a(test_data_dir() + "oco2_sim_scene.h5", sid);
  BOOST_CHECK_CLOSE(a.co2_vmr(77731.75069322), 0.000377794, 1e-4);
  BOOST_CHECK_CLOSE(a.co2_vmr(86876.6640625), 0.00037741665, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
