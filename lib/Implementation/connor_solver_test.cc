#include "connor_solver.h"
#include "chisq_convergence.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(connor_solver, GlobalFixture)

// This is a simple nonlinear problem that has the solution 1, 2, 3,
// 4, 5.
class CostFunctionTest : public CostFunction {
public:
  virtual void cost_function(const blitz::Array<double, 1>& x,
			blitz::Array<double, 1>& Residual,
			blitz::Array<double, 1>& Se,
			blitz::Array<double, 2>& Jacobian) const
  {
    Jacobian.resize(10, 5);
    Jacobian = 0;
    Jacobian(0, 0) = 2 * x(0);
    Jacobian(1, 1) = 2 * x(1);
    Jacobian(2, 2) = 2 * x(2);
    Jacobian(3, 3) = 2 * x(3);
    Jacobian(4, 4) = 2 * x(4);
    Jacobian(5,0) = 1.0;
    Jacobian(5,1) = 1.0;
    Jacobian(6,1) = 1.0;
    Jacobian(6,2) = 1.0;
    Jacobian(7,1) = 2.0;
    Jacobian(7,3) = 1.0;
    Jacobian(8,3) = 2.0;
    Jacobian(8,4) = -1.0;
    Jacobian(9,2) = 2.0;
    Jacobian(9,3) = 1.0;
    Residual.resize(10);
    Residual = (x(0) * x(0) - 1), 
      (x(1) * x(1) - 4),
      (x(2) * x(2) - 9),
      (x(3) * x(3) - 16),
      (x(4) * x(4) - 25),
      (x(0) + x(1) - 3),
      (x(1) + x(2) - 5),
      (2 * x(1) + x(3) - 8),
      (2 * x(3) - x(4) - 3),
      (2 * x(2) + x(3) - 10);
    Se.resize(10);
    Se = 1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6;
  }
  virtual void print(std::ostream& Os) const
  { Os << "CostFunctionTest"; }
};

// This is CostFunctionTest but with an extra ignored parameter. This
// should go to the apriori and stick there.
class CostFunctionRankDegenerate : public CostFunction {
public:
  virtual void cost_function(const blitz::Array<double, 1>& x,
			blitz::Array<double, 1>& Residual,
			blitz::Array<double, 1>& Se,
			blitz::Array<double, 2>& Jacobian) const
  {
    Jacobian.resize(10, 6);
    Jacobian = 0;
    Jacobian(0, 0) = 2 * x(0);
    Jacobian(1, 1) = 2 * x(1);
    Jacobian(2, 2) = 2 * x(2);
    Jacobian(3, 3) = 2 * x(3);
    Jacobian(4, 4) = 2 * x(4);
    Jacobian(5,0) = 1.0;
    Jacobian(5,1) = 1.0;
    Jacobian(6,1) = 1.0;
    Jacobian(6,2) = 1.0;
    Jacobian(7,1) = 2.0;
    Jacobian(7,3) = 1.0;
    Jacobian(8,3) = 2.0;
    Jacobian(8,4) = -1.0;
    Jacobian(9,2) = 2.0;
    Jacobian(9,3) = 1.0;
    Residual.resize(10);
    Residual = (x(0) * x(0) - 1), 
      (x(1) * x(1) - 4),
      (x(2) * x(2) - 9),
      (x(3) * x(3) - 16),
      (x(4) * x(4) - 25),
      (x(0) + x(1) - 3),
      (x(1) + x(2) - 5),
      (2 * x(1) + x(3) - 8),
      (2 * x(3) - x(4) - 3),
      (2 * x(2) + x(3) - 10);
    Se.resize(10);
    Se = 1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6;
  }
  virtual void print(std::ostream& Os) const
  { Os << "CostFunctionTest"; }
};

BOOST_AUTO_TEST_CASE(basic)
{
  ConnorSolver cs(boost::shared_ptr<CostFunction>(new CostFunctionTest),
		  boost::shared_ptr<ConvergenceCheck>(new ChisqConvergence));
  Array<double, 1> initial_guess(5);
  Array<double, 1> apriori(5);
  Array<double, 2> apriori_cov(5, 5);
  apriori_cov = 0;
  apriori_cov(0,0) = 1;
  apriori_cov(1,1) = 1;
  apriori_cov(2,2) = 1;
  apriori_cov(3,3) = 1;
  apriori_cov(4,4) = 1;
  apriori = 1,2,3,4,5;
  initial_guess = 2, 3, 4, 5, 6;
  cs.solve(initial_guess, apriori, apriori_cov);
  BOOST_CHECK_CLOSE(cs.x_solution()(0), 1.0, 1e-2);
  BOOST_CHECK_CLOSE(cs.x_solution()(1), 2.0, 1e-2);
  BOOST_CHECK_CLOSE(cs.x_solution()(2), 3.0, 1e-2);
  BOOST_CHECK_CLOSE(cs.x_solution()(3), 4.0, 1e-2);
  BOOST_CHECK_CLOSE(cs.x_solution()(4), 5.0, 1e-2);
}

BOOST_AUTO_TEST_CASE(rank_deficient)
{
  ConnorSolver cs(boost::shared_ptr<CostFunction>(new 
						  CostFunctionRankDegenerate),
		  boost::shared_ptr<ConvergenceCheck>(new ChisqConvergence));
  Array<double, 1> initial_guess(6);
  Array<double, 1> apriori(6);
  Array<double, 2> apriori_cov(6, 6);
  apriori_cov = 0;
  apriori_cov(0,0) = 1;
  apriori_cov(1,1) = 1;
  apriori_cov(2,2) = 1;
  apriori_cov(3,3) = 1;
  apriori_cov(4,4) = 1;
  apriori_cov(5,5) = 1e-6;
  apriori = 1,2,3,4,5,6;
  initial_guess = 2, 3, 4, 5, 6, 6;
  cs.solve(initial_guess, apriori, apriori_cov);
  BOOST_CHECK_CLOSE(cs.x_solution()(0), 1.0, 1e-2);
  BOOST_CHECK_CLOSE(cs.x_solution()(1), 2.0, 1e-2);
  BOOST_CHECK_CLOSE(cs.x_solution()(2), 3.0, 1e-2);
  BOOST_CHECK_CLOSE(cs.x_solution()(3), 4.0, 1e-2);
  BOOST_CHECK_CLOSE(cs.x_solution()(4), 5.0, 1e-2);
  BOOST_CHECK_CLOSE(cs.x_solution()(5), 6.0, 1e-2);
}

BOOST_AUTO_TEST_CASE(do_inversion_test)
{
  is_long_test();		// Skip unless we are running long tests.

  // We get the expected results here by running the old Fortran
  // version and capturing the calculated values. We make sure that we
  // are calculating the same values with the new code.
  ConnorSolver cs(boost::shared_ptr<CostFunction>(new 
						  CostFunctionRankDegenerate),
		  boost::shared_ptr<ConvergenceCheck>(new ChisqConvergence));
  blitz::Array<double, 1> dx;
  blitz::Array<double, 2> kt_se_m1_k;
  cs.test_do_inversion(test_data_dir() + 
		       "expected/connor_solver/connor_save.txt", dx, 
		       kt_se_m1_k);
  FitStatistic fstat = cs.fit_statistic();
  BOOST_CHECK_CLOSE(fstat.d_sigma_sq_scaled, 2898.0253091504742, 1e-4);
  BOOST_CHECK_CLOSE(fstat.d_sigma_sq, 310088.73232535017, 1e-4);
  BOOST_CHECK_CLOSE(fstat.gamma2(), 310323.31238974462, 1e-4);
  BOOST_CHECK_CLOSE(fstat.gamma2_fc(), 68.272218668209760, 1e-4);
  BOOST_CHECK_CLOSE(fstat.chisq_apriori, 0.0, 1e-4);
  BOOST_CHECK_CLOSE(fstat.chisq_measured, 310323.31238974462, 1e-4);
  BOOST_CHECK_CLOSE(fstat.chisq_apriori_fc, 16.630739518037689, 1e-4);
  BOOST_CHECK_CLOSE(fstat.chisq_measured_fc, 51.641434095666973, 1e-4);
  BOOST_CHECK_CLOSE(dx(10), -1.9295747636982488e-06, 1e-4);
  BOOST_CHECK_CLOSE(dx(20), -493.30933006057325, 1e-4);
  BOOST_CHECK_CLOSE(dx(21), -4.5497933230390259, 1e-4);
  BOOST_CHECK_CLOSE(dx(30), 0.63245395057915621, 1e-4);
  BOOST_CHECK_CLOSE(kt_se_m1_k(0,0), 29584028.051484555, 1e-4);
  BOOST_CHECK_CLOSE(kt_se_m1_k(0,1), 55199185.339754276, 1e-4);
  BOOST_CHECK_CLOSE(kt_se_m1_k(1,1), 108736851.88356204, 1e-4);
}

BOOST_AUTO_TEST_CASE(averaging_kernel_test)
{
  // This data comes from running a test case to the end, and then
  // saving the state. We restore the state to ConnorSolver (so it is
  // like we just finished this long convergence calculation), and
  // check that we then calculate the averaging kernel correctly.
  ConnorSolver cs(boost::shared_ptr<CostFunction>(new 
						  CostFunctionRankDegenerate),
		  boost::shared_ptr<ConvergenceCheck>(new ChisqConvergence));
  IfstreamCs in(test_data_dir() + "connor_converged.txt");
  in >> cs;
  Array<double, 2> ak(cs.averaging_kernel());
  BOOST_CHECK_CLOSE(ak(0,0), 2.0114592943839052e-05, 2e-3);
  BOOST_CHECK_CLOSE(ak(0,1), 2.7111012604872633e-05, 2e-3);
  BOOST_CHECK_CLOSE(ak(1,1), 0.00044724517297504583, 2e-3);
  BOOST_CHECK_CLOSE(ak(1,2), 0.00069114989773355032, 2e-3);
  BOOST_CHECK_CLOSE(ak(106,106), 0.99991041825545657, 2e-3);
}

// This data comes from a test with a nondiagonal covariance
// matrix. This checks the correct handling of the apriori 
// covariance matrix with a degenerate jacobian
BOOST_AUTO_TEST_CASE(nondiagonal_cov)
{
  is_long_test();		// Skip unless we are running long tests.

  // We get the expected results here by running the old Fortran
  // version and capturing the calculated values. We make sure that we
  // are calculating the same values with the new code.
  ConnorSolver cs(boost::shared_ptr<CostFunction>(new 
						  CostFunctionRankDegenerate),
		  boost::shared_ptr<ConvergenceCheck>(new ChisqConvergence));
  blitz::Array<double, 1> dx;
  blitz::Array<double, 2> kt_se_m1_k;
  cs.test_do_inversion(test_data_dir() + 
		       "expected/connor_solver/connor_nondiagonal_cov.txt", dx, 
		       kt_se_m1_k);
  BOOST_CHECK_CLOSE(dx(10), 1.80693604E-05, 1e-2);
}

BOOST_AUTO_TEST_SUITE_END()
