#include "unit_test_support.h"
#include "fts_run_log_output.h"
#include "output_hdf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(fts_run_log_output, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  FtsRunLog rl(test_data_dir() + "in/tccon_runlog.grl");
  std::vector<FtsRunLogRecord> rlr_vec;
  rlr_vec.push_back(rl.read("pa20091103saaaaa_100223160344.008"));

  FtsRunLogOutput po(rlr_vec);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("fts_run_log_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("fts_run_log_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


