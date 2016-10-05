#include "unit_test_support.h"
#include "fp_exception.h"
#include "solver_finished_fixture.h"
#include "output_hdf.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(full_output_test, SolverFinishedFixture)

BOOST_AUTO_TEST_CASE(full_output)
{
  is_long_test();		// Skip unless we are running long tests.

  boost::shared_ptr<OutputHdf> out(new OutputHdf("full_output.h5", 20, 112, 5, 3));
//  add_file_to_cleanup("configuration_heritage_output.h5");
  config_state_vector->update_state(initial_sv, apriori_cov);
  BOOST_FOREACH(boost::shared_ptr<RegisterOutputBase> i, 
		config_register_output) {
    i->register_output(out);
    i->register_output_apriori(out);
  }
  config_state_vector->update_state(solver->x_solution(), 
				    solver->aposteriori_covariance());
  out->write();
}  
BOOST_AUTO_TEST_SUITE_END()
