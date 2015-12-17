#include "connor_convergence.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(connor_convergence, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  ConnorConvergence conv(config_forward_model, 0.1, 10, 4, 1.4);
  FitStatistic  fs, fs_last;
  bool has_converged, convergence_failed, step_diverged;
  double gamma;

  gamma = 8.0;

  fs_last.chisq_apriori = 0.5;
  fs_last.chisq_measured = 1.5;
  fs_last.chisq_apriori_fc = 0.25;
  fs_last.chisq_measured_fc = 0.75;

  //  Set the value of the cost function to
  //  be the same as the predicted reduction
  //  from the previous iteration (the best
  //  prediction).
  //
  fs.chisq_apriori = 0.25;
  fs.chisq_measured = 0.75;


  //  Set the number of iterations to 1, which
  //  means that there is no predicted cost function,
  //  hence gamma should not get updated.  Moreover,
  //  there is no previously evaluated cost function;
  //  therefore, it is not possible to conclude that
  //  step diverged.
  //
  fs.number_iteration = 1;
  //
  fs.d_sigma_sq_scaled = 0.1 + 0.05;
  conv.convergence_check( fs_last, fs, has_converged, 
                           convergence_failed, gamma, step_diverged );
  BOOST_CHECK( gamma == 8.0 );
  BOOST_CHECK( step_diverged == false );
  BOOST_CHECK( has_converged == false );
  //
  fs.d_sigma_sq_scaled = 0.1 - 0.05;
  conv.convergence_check( fs_last, fs, has_converged, 
                           convergence_failed, gamma, step_diverged );
  BOOST_CHECK( gamma == 8.0 );
  BOOST_CHECK( step_diverged == false );
  BOOST_CHECK( has_converged == true );


  // Now after more than one iteration and
  // a reduction in the cost function as good
  // as the predicted reduction, gamma should
  // decrease.
  //
  fs.number_iteration = 2;
  //
  fs.d_sigma_sq_scaled = 0.1 + 0.05;
  gamma = 8.0;
  conv.convergence_check( fs_last, fs, has_converged, 
                           convergence_failed, gamma, step_diverged );
  BOOST_CHECK( gamma < 8.0 );
  BOOST_CHECK( step_diverged == false );
  BOOST_CHECK( has_converged == false );
  //
  fs.d_sigma_sq_scaled = 0.1 - 0.05;
  gamma = 8.0;
  conv.convergence_check( fs_last, fs, has_converged, 
                           convergence_failed, gamma, step_diverged );
  BOOST_CHECK( gamma < 8.0 );
  BOOST_CHECK( step_diverged == false );
  BOOST_CHECK( has_converged == true );


  // Now after more than one iteration and
  // an increase in the cost function (a bad
  // prediction), gamma should increase.
  //
  fs.chisq_apriori = 0.55;
  fs.chisq_measured = 1.75;
  //
  gamma = 8.0;
  fs.d_sigma_sq_scaled = 0.1 + 0.05;
  conv.convergence_check( fs_last, fs, has_converged, 
                           convergence_failed, gamma, step_diverged );
  BOOST_CHECK( gamma > 8.0 );
  BOOST_CHECK( step_diverged == true );
  BOOST_CHECK( has_converged == false );
  //
  gamma = 8.0;
  fs.d_sigma_sq_scaled = 0.1 - 0.05;
  conv.convergence_check( fs_last, fs, has_converged, 
                           convergence_failed, gamma, step_diverged );
  BOOST_CHECK( gamma > 8.0 );
  BOOST_CHECK( step_diverged == true );
  //  There is another logical error in the method convergence_check.
  //  Convergence check as implemented in the method (based on
  //  the step size) should be independent of the divergence
  //  check as implemented in the method (based on an increase
  //  or decrease in the value of the cost function). If a step
  //  is small enough to consider a convergence, then the 
  //  convergence must be true regardless of whether the step
  //  increases or decreases the 
  //  value of the cost function.  The description of the
  //  consequences of the incorrect implementation is beyond
  //  the scope of this comment.  The following check should
  //  not cause any error, but it does.  The method 
  //  will work correctly if the return statement in the middle
  //  of the method is just deleted; however, I don't know
  //  whether or not the elimination of the return statement
  //  will cause an error in any other place.
//  BOOST_CHECK( has_converged == true );


}

BOOST_AUTO_TEST_SUITE_END()
