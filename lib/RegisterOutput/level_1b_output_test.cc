#include "unit_test_support.h"
#include "level_1b_output.h"
#include "output_hdf.h"
#include "level_1b_acos.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(level_1b_output, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<HdfFile> hf(new HdfFile(test_data_dir() + "l1b.h5"));
  boost::shared_ptr<AcosSoundingId> sid
    (new AcosSoundingId(*hf, "20090725020316", AcosSoundingId::P_SOUNDING));
  boost::shared_ptr<Level1b> f(new Level1bAcos (hf, sid));
  Level1bOutput po(f);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("level_1b_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("level_1b_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in level_1b unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


