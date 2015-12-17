#include "unit_test_support.h"
#include "connor_solver_output.h"
#include "output_hdf.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(connor_solver_output, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  // This data comes from running a test case to the end, and then
  // saving the state. We restore the state to ConnorSolver (so it is
  // like we just finished this long convergence calculation), and
  // check that we then calculate the averaging kernel correctly.
  boost::shared_ptr<CostFunction> cf;
  boost::shared_ptr<ConvergenceCheck> ccheck;
  boost::shared_ptr<ConnorSolver> cs(new ConnorSolver(cf, ccheck));
  IfstreamCs in(test_data_dir() + "connor_converged.txt");
  in >> *cs;
  ConnorSolverOutput conout(cs);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("connor_solver_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("connor_solver_output.h5");
  conout.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculations are checked in ConnorSolver unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()



