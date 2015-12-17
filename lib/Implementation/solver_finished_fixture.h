#ifndef SOLVER_FINISHED_FIXTURE_H
#define SOLVER_FINISHED_FIXTURE_H
#include "global_fixture.h"
#include "configuration_fixture.h"

namespace FullPhysics {
/****************************************************************//**
  This is a test fixture that puts us into the state we are in after
  doing the solver step of l2_fp. This is for the same run as
  l2_fp_run does, we just skip over the actual solving step. This is
  used to test those things that occur at the end of the l2_fp
  process, e.g., the output.

  This uses the file "connor_converged.txt" found in unit_test_data in
  the source directory. You may need to periodically regenerate this
  data as we change the Level 2 full physics algorithms. This can be
  generated using l2_fp with the "-t" option (see l2_fp -h for
  details). 
*******************************************************************/
class SolverFinishedFixture : public ConfigurationFixture {
public:
  SolverFinishedFixture();
  virtual ~SolverFinishedFixture() {}
  boost::shared_ptr<ConnorSolver> solver;
  blitz::Array<double, 1> initial_sv;
  blitz::Array<double, 1> apriori_sv;
  blitz::Array<double, 2> apriori_cov;
};
}
#endif
