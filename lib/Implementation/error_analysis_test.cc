#include "error_analysis.h"
#include "absorber_absco.h"
#include "unit_test_support.h"
#include "solver_finished_fixture.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(error_analysis, SolverFinishedFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  is_long_test();		// Skip unless we are running long tests.
  IfstreamCs expected_data(test_data_dir() + 
			   "expected/error_analysis/basic");
  const ErrorAnalysis& err = *config_error_analysis;
  BOOST_CHECK_CLOSE(err.signal(0), 3.6830670117200004e-07, 1e-2);
  BOOST_CHECK_CLOSE(err.signal(1), 7.6862606524800006e-08, 1e-2);
  BOOST_CHECK_CLOSE(err.signal(2), 2.6669225391799998e-08, 1e-2);
  BOOST_CHECK_CLOSE(err.noise(0), 1.9218877132799996e-09, 1e-2);
  BOOST_CHECK_CLOSE(err.noise(1), 1.26415597811e-09, 1e-2);
  BOOST_CHECK_CLOSE(err.noise(2), 5.939573830669999e-10, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_sum_sq(0), 2.4365234387594381e-16, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_sum_sq(1), 3.3238014086328966e-18,
		    1e-2);
  BOOST_CHECK_CLOSE(err.residual_sum_sq(2), 8.7843864248064908e-18,
		    1e-2);
  BOOST_CHECK_CLOSE(err.residual_mean_sq(0), 4.5004141665502041e-10, 1e-2);
  BOOST_CHECK_CLOSE(err.residual_mean_sq(1), 7.4367006091835221e-11, 
		    1e-2);
  BOOST_CHECK_CLOSE(err.residual_mean_sq(2), 1.4194250697466709e-10, 
		    1e-2);
  BOOST_CHECK_CLOSE(err.relative_residual_mean_sq(0), 
		    err.residual_mean_sq(0) / err.signal(0), 1e-2);
  BOOST_CHECK_CLOSE(err.relative_residual_mean_sq(1), 
		    err.residual_mean_sq(1) / err.signal(1), 
		    1e-2);
  BOOST_CHECK_CLOSE(err.relative_residual_mean_sq(2), 
	    err.residual_mean_sq(2) / err.signal(2), 
	    1e-2);
  BOOST_CHECK_CLOSE(err.reduced_chisq(0), 0.05402455082526629, 1e-2);
  BOOST_CHECK_CLOSE(err.reduced_chisq(2), 0.054857252629198799,
		    1e-2);
  BOOST_CHECK_CLOSE(err.reduced_chisq(1), 0.0034802312226315355,
		    1e-2);
  BOOST_CHECK_CLOSE(err.xco2_uncert_noise(), 
		    sqrt(err.xco2_measurement_error()), 1e-4);
  BOOST_CHECK_CLOSE(err.xco2_uncert_smooth(), 
		    sqrt(err.xco2_smoothing_error()), 1e-4);
  BOOST_CHECK_CLOSE(err.xco2_uncert_interf(), 
		    sqrt(err.xco2_interference_error()), 1e-4);
  BOOST_CHECK_CLOSE(err.xco2_measurement_error(), 8.0542354367629385e-12,
		    1e-4);
  BOOST_CHECK_CLOSE(err.xco2_smoothing_error(), 1.6529451260379834e-12, 1e-4);
  BOOST_CHECK_CLOSE(err.xco2_interference_error(), 4.1718475757889482e-12,
		    1e-4);
  // Equation 3-109 in the ATB
  BOOST_CHECK_CLOSE(sqr(err.xco2_uncertainty()),
		    err.xco2_measurement_error() +
		    err.xco2_smoothing_error() +
		    err.xco2_interference_error(), 1e-4);
  BOOST_CHECK_CLOSE(err.degrees_of_freedom_full_vector(), 
		    15.762273734062212, 1e-4);
  BOOST_CHECK_CLOSE(err.degrees_of_freedom_xco2(), 
		    0.91615188124993685, 1e-4);

  Array<double, 1> xco2_avg_kernel_expected, xco2_avg_kernel_norm_expected,
    xco2_correlation_interf_expected, isu_expected, xco2_gain_expected;
  expected_data >> xco2_avg_kernel_expected
		>> xco2_avg_kernel_norm_expected
		>> xco2_correlation_interf_expected
		>> isu_expected
        >> xco2_gain_expected;
  BOOST_CHECK_MATRIX_CLOSE(err.xco2_avg_kernel(), xco2_avg_kernel_expected);
  BOOST_CHECK_MATRIX_CLOSE_TOL(err.xco2_avg_kernel_norm(), 
			       xco2_avg_kernel_norm_expected, 1e-6);
  BOOST_CHECK_MATRIX_CLOSE(err.xco2_correlation_interf(), 
			   xco2_correlation_interf_expected);
  BOOST_CHECK_MATRIX_CLOSE_TOL(err.interference_smoothing_uncertainty(),
			       isu_expected, 1e-13);
  BOOST_CHECK_MATRIX_CLOSE_TOL(err.xco2_gain_vector(), xco2_gain_expected, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
