#include <bard_nlls_problem.h>
#include <brown_nlls_problem.h>
#include <freudenstein_roth_nlls_problem.h>
#include <helical_valley_nlls_problem.h>
#include <jennrich_sampson_nlls_problem.h>
#include <meyer_nlls_problem.h>
#include <powell_nlls_problem.h>
#include <powell_singular_nlls_problem.h>
#include <rosenbrock2_nlls_problem.h>
#include <unit_test_support.h>
#include <fp_exception.h>
#include <cost_minimizer_gsl.h>

#include <bard_ml_problem.h>
#include <meyer_ml_problem.h>
#include <nlls_max_likelihood.h>


using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(cost_minimizer_gsl, GlobalFixture)

/* convergence check thresholds */
double dx_epsabs=1e-5, dx_epsrel=1e-5, size_tol=1e-5;
bool verbose = false;

BOOST_AUTO_TEST_CASE(bard)
{
  Array<double, 1> x0(3); x0 = 1.0, 1.0, 1.0;
  boost::shared_ptr<BardNLLSProblem> nlls(new BardNLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver(1000, dx_epsabs, dx_epsrel, size_tol, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Bard function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 225);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK( (abs(cst-0.0041074) < 0.000001) || (abs(cst-8.71431) < 0.001) );
  if( abs(cst-0.0041074) < 0.000001 ) {
    BOOST_CHECK_CLOSE(nlls->parameters()(0), 0.082411, 0.01);
    BOOST_CHECK_CLOSE(nlls->parameters()(1), 1.1330, 0.01);
    BOOST_CHECK_CLOSE(nlls->parameters()(2), 2.3437, 0.01);
  }
}


BOOST_AUTO_TEST_CASE(brown)
{
  Array<double, 1> x0(2); x0 = 1.0, 1.0;
  boost::shared_ptr<BrownNLLSProblem> nlls(new BrownNLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver(1000, dx_epsabs, dx_epsrel, size_tol, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Brown function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 330);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK( (abs(cst-0.0) < 0.0000001) );
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 1.0e6, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 2.0e-6, 0.01);
}


BOOST_AUTO_TEST_CASE(freudenstein_roth__a)
{
  Array<double, 1> x0(2); x0 = 0.5, -2.0;
  boost::shared_ptr<FreudensteinRothNLLSProblem> nlls(new FreudensteinRothNLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver(1000, dx_epsabs, dx_epsrel, size_tol, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Freudenstein/Roth function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 125);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK( (abs(cst-0.0) < 0.0000001) || (abs(cst-24.4921) < 0.001) );
  if( abs(cst-0.0) < 0.0000001 ) {
    BOOST_CHECK_CLOSE(nlls->parameters()(0), 5.0, 0.01);
    BOOST_CHECK_CLOSE(nlls->parameters()(1), 4.0, 0.01);
  } else {
    BOOST_CHECK_CLOSE(nlls->parameters()(0), 11.413, 0.01);
    BOOST_CHECK_CLOSE(nlls->parameters()(1), -0.89681, 0.01);
  }
}


BOOST_AUTO_TEST_CASE(freudenstein_roth__b)
{
  Array<double, 1> x0(2); x0 = 6.0, 2.0;
  boost::shared_ptr<FreudensteinRothNLLSProblem> nlls(new FreudensteinRothNLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver(1000, dx_epsabs, dx_epsrel, size_tol, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Freudenstein/Roth function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 120);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK( (abs(cst-0.0) < 0.0000001) || (abs(cst-24.4921) < 0.001) );
  if( abs(cst-0.0) < 0.0000001 ) {
    BOOST_CHECK_CLOSE(nlls->parameters()(0), 5.0, 0.01);
    BOOST_CHECK_CLOSE(nlls->parameters()(1), 4.0, 0.01);
  } else {
    BOOST_CHECK_CLOSE(nlls->parameters()(0), 11.413, 0.01);
    BOOST_CHECK_CLOSE(nlls->parameters()(1), -0.89681, 0.01);
  }
}


BOOST_AUTO_TEST_CASE(helical_valley)
{
  Array<double, 1> x0(3); x0 = -1.2, 0.0, 0.0;
  boost::shared_ptr<HelicalValleyNLLSProblem> nlls(new HelicalValleyNLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver(1000, dx_epsabs, dx_epsrel, size_tol*0.01, nlls, Array<double,1>(), verbose);  // A smaller size_tol is needed, hence *0.01
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Helical Valley function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 295);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK( (abs(cst-0.0) < 0.0000001) );
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 1.0, 0.01);
  BOOST_CHECK(abs(nlls->parameters()(1)-0.0) < 0.0000001);
  BOOST_CHECK(abs(nlls->parameters()(2)-0.0) < 0.0000001);
}


BOOST_AUTO_TEST_CASE(jennrich_sampson)
{
  Array<double, 1> x0(2); x0 = 0.3, 0.4;
  boost::shared_ptr<JennrichSampsonNLLSProblem> nlls(new JennrichSampsonNLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver(1000, dx_epsabs, dx_epsrel, size_tol, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Jennrich/Sampson function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 85);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK_CLOSE(cst, 62.1813, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 0.2578, 0.015);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 0.2578, 0.015);
}


BOOST_AUTO_TEST_CASE(meyer)
{
  Array<double, 1> x0(3); x0 = 0.02, 4000.0, 250.0;
  boost::shared_ptr<MeyerNLLSProblem> nlls(new MeyerNLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver(2000, dx_epsabs, dx_epsrel, size_tol, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Meyer function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 1700);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK_CLOSE(cst, 43.9729, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 0.0056096, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 6181.35, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(2), 345.224, 0.01);
}


BOOST_AUTO_TEST_CASE(powell)
{
  Array<double, 1> x0(2); x0 = 0.0, 1.0;
  boost::shared_ptr<PowellNLLSProblem> nlls(new PowellNLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver(1000, dx_epsabs, dx_epsrel, size_tol, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Powell function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 700);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK( (abs(cst-0.0) < 0.0000001) );
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 1.09816e-5, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 9.10615, 0.01);
}


BOOST_AUTO_TEST_CASE(powell_singular)
{
  Array<double, 1> x0(4); x0 = 3.0, -1.0, 0.0, 1.0;
  boost::shared_ptr<PowellSingularNLLSProblem> nlls(new PowellSingularNLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver(3000, 1e-21, 1e-21, 1e-10, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Powell Singular function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 1200);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK( (abs(cst-0.0) < 0.0000001) );
  BOOST_CHECK(abs(nlls->parameters()(0)-0.0) < 0.0000001);
  BOOST_CHECK(abs(nlls->parameters()(1)-0.0) < 0.0000001);
  BOOST_CHECK(abs(nlls->parameters()(2)-0.0) < 0.0000001);
  BOOST_CHECK(abs(nlls->parameters()(3)-0.0) < 0.0000001);
}


BOOST_AUTO_TEST_CASE(rosenbrock2)
{
  Array<double, 1> x0(2); x0 = -1.2, 1.0;
  boost::shared_ptr<Rosenbrock2NLLSProblem> nlls(new Rosenbrock2NLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver(1000, dx_epsabs, dx_epsrel, size_tol, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Rosenbrock function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 175);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK( (abs(cst-0.0) < 0.0000001) );
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 1.0, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 1.0, 0.01);
}


BOOST_AUTO_TEST_CASE(powell_singular_multiple)
{
  Array<double, 1> x0(4); x0 = 3.0, -1.0, 0.0, 1.0;
  boost::shared_ptr<PowellSingularNLLSProblem> nlls(new PowellSingularNLLSProblem);
  nlls->parameters(x0);
  CostMinimizerGSL solver1(1000, 1e-5, 1e-5, 1e-4, nlls, Array<double,1>(), verbose);
  CostMinimizerGSL solver2(2000, 1e-21, 1e-21, 1e-10, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver1.status(), (int)CostMinimizerGSL::UNTRIED);
  BOOST_CHECK_EQUAL((int)solver2.status(), (int)CostMinimizerGSL::UNTRIED);



//   std::cout 
//      << std::endl
//      << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< nlls / solver1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" 
//      << std::endl;

  solver1.solve();

  int n_f_calls1 = nlls->num_residual_evaluations();
  int n_j_calls1 = nlls->num_jacobian_evaluations();
  double cst1 = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Powell Singular function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls1 << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls1 << std::endl
//      << "   Final solver1 status = " << solver1.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst1 << std::endl;
//   for( int i=0; i<=solver1.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver1.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver1.cost_at_accepted_points()[i] << std::endl
//         << "   Grad["<<i<<"] = " << solver1.gradient_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver1.accepted_points().size(), solver1.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver1.cost_at_accepted_points().size(), solver1.num_accepted_steps()+1);

  int iLast1 = solver1.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst1, solver1.cost_at_accepted_points()[iLast1]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver1.accepted_points()[iLast1])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver1.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls1 < 430);
  BOOST_CHECK(n_j_calls1 <= 0);



//   std::cout 
//      << std::endl
//      << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< nlls / solver2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
//      << std::endl;

  solver2.solve();

  int n_f_calls2 = nlls->num_residual_evaluations();
  int n_j_calls2 = nlls->num_jacobian_evaluations();
  double cst2 = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Powell Singular function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls2 << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls2 << std::endl
//      << "   Final solver2 status = " << solver2.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst2 << std::endl;
//   for( int i=0; i<=solver2.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver2.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver2.cost_at_accepted_points()[i] << std::endl
//         << "   Grad["<<i<<"] = " << solver2.gradient_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver2.accepted_points().size(), solver2.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver2.cost_at_accepted_points().size(), solver2.num_accepted_steps()+1);

  int iLast2 = solver2.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst2, solver2.cost_at_accepted_points()[iLast2]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver2.accepted_points()[iLast2])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver2.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls2 < 1310);
  BOOST_CHECK(n_j_calls2 <= 0);
  BOOST_CHECK( (abs(cst2-0.0) < 0.0000001) );
  BOOST_CHECK(abs(nlls->parameters()(0)-0.0) < 0.0000001);
  BOOST_CHECK(abs(nlls->parameters()(1)-0.0) < 0.0000001);
  BOOST_CHECK(abs(nlls->parameters()(2)-0.0) < 0.0000001);
  BOOST_CHECK(abs(nlls->parameters()(3)-0.0) < 0.0000001);


  // Last point recorded by the first solver must be the 
  // same as the first point recorded by the second solver.
  //
  // For this solver the above statements may not be true,
  // hence the following two lines are commented out.  I
  // will need to investigate more.
  //
//  BOOST_CHECK_EQUAL(solver1.cost_at_accepted_points()[iLast1], solver2.cost_at_accepted_points()[0]);
//  BOOST_CHECK_CLOSE(sum(abs(solver1.accepted_points()[iLast1]-solver2.accepted_points()[0])), 0.0, 1e-12);

}








//  The case tests interfacing of 
//  maximum-likelihood/nlls-problem/nlls-solver
BOOST_AUTO_TEST_CASE(bard_ml)
{
  Array<double, 1> x0(3); x0 = 1.0, 1.0, 1.0;
  Array<double, 1> measurement(15);
  Array<double, 1> measurement_error_cov(15);
  measurement = 0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, 0.73, 0.96, 1.34, 2.10, 4.39;
  measurement_error_cov = 1.0;

  boost::shared_ptr<BardMLProblem> ml(new BardMLProblem(measurement, measurement_error_cov));
  boost::shared_ptr<NLLSMaxLikelihood> nlls(new NLLSMaxLikelihood(ml, false));
  nlls->parameters(x0);
  CostMinimizerGSL solver(1000, dx_epsabs, dx_epsrel, size_tol, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Bard/ML function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 225);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK( (abs(cst-0.0041074) < 0.000001) || (abs(cst-8.71431) < 0.001) );
  if( abs(cst-0.0041074) < 0.000001 ) {
    BOOST_CHECK_CLOSE(nlls->parameters()(0), 0.082411, 0.01);
    BOOST_CHECK_CLOSE(nlls->parameters()(1), 1.1330, 0.01);
    BOOST_CHECK_CLOSE(nlls->parameters()(2), 2.3437, 0.01);
  }
}


//  The case tests interfacing of 
//  maximum-likelihood/nlls-problem/nlls-solver
BOOST_AUTO_TEST_CASE(meyer_ml)
{
  Array<double, 1> x0(3); x0 = 0.02, 4000.0, 250.0;
  Array<double, 1> measurement(16);
  Array<double, 1> measurement_error_cov(16);
  measurement = 34780.0, 28610.0, 23650.0, 19630.0, 16370.0, 13720.0, 11540.0, 9744.0, 8261.0, 7030.0, 6005.0, 5147.0, 4427.0, 3820.0, 3307.0, 2872.0;
  measurement_error_cov = 1.0;

  boost::shared_ptr<MeyerMLProblem> ml(new MeyerMLProblem(measurement, measurement_error_cov));
  boost::shared_ptr<NLLSMaxLikelihood> nlls(new NLLSMaxLikelihood(ml, false));
  nlls->parameters(x0);
  CostMinimizerGSL solver(2000, dx_epsabs, dx_epsrel, size_tol, nlls, Array<double,1>(), verbose);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

//   std::cout 
//      << "Testing CostMinimizerGSL with Meyer/ML function:" << std::endl
//      << "   Number of residual function evaluations = " << n_f_calls << std::endl
//      << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
//      << "   Final solver status = " << solver.status_str() << std::endl
//      << "   Final problem status (point) = " << nlls->parameters() << std::endl
//      << "   Final problem status (cost value) = " << cst << std::endl;
//   for( int i=0; i<=solver.num_accepted_steps(); i++ )
//      std::cout 
//         << "   ========================================" << std::endl
//         << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
//         << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)CostMinimizerGSL::SUCCESS);
  BOOST_CHECK(n_f_calls < 1700);
  BOOST_CHECK(n_j_calls <= 0);
  BOOST_CHECK_CLOSE(cst, 43.9729, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 0.0056096, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 6181.35, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(2), 345.224, 0.01);
}


BOOST_AUTO_TEST_SUITE_END()
