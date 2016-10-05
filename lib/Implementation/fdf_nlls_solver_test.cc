#include <gsl/gsl_blas.h>
#include <gsl_lsp.h>
#include <fdf_nlls_solver.h>
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



using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(gsl_fdf_nlls_solver, GlobalFixture)


/* for selecting a problem */
gsl_multifit_function_fdf f;

/* for initial guess */
gsl_vector *x0 = (gsl_vector *)0;

/* convergence check thresholds */
double dx_epsabs=1e-5, dx_epsrel=1e-5, g_epsabs=1e-5;

/* select the solver */
const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmsder;
/*const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmder;*/

/* for final result */
gsl_vector *x = (gsl_vector *)0;

/* to optionally get residual at the final point */
gsl_vector *r = (gsl_vector *)0;

/* to optionally get Jacobian at the final point */
gsl_matrix *J = (gsl_matrix *)0;

/* to get the number of iterations */
unsigned int num_iter;

/* to get the GSL solver status */
int slvr_status;

/*  */
short int verbose = false;


BOOST_AUTO_TEST_CASE(bard)
{
    BardNLLSProblem bard_nlls;

    /* selection of the problem */
    f = gsl_get_lsp_fdf(&bard_nlls);

    /* for the initial guess */
    x0 = gsl_vector_calloc(f.p);

    /* comment out any of the following 3 lines if the optional 
     * solution, residual and Jacobian at the final point are not
     * desired. */
    x = gsl_vector_alloc(f.p);
    r = gsl_vector_alloc(f.n);
    /*J = gsl_matrix_alloc(f.n,f.p);*/

    /* selection of the initial guess */
    gsl_vector_set(x0, 0, 1.0);
    gsl_vector_set(x0, 1, 1.0);
    gsl_vector_set(x0, 2, 1.0);

    /* solve */
    slvr_status = fdf_nlls_solver(&f, x0, T, dx_epsabs, dx_epsrel, g_epsabs, 100, verbose, &num_iter, x, r, J);
    int n_f_calls = bard_nlls.num_residual_evaluations();
    int n_j_calls = bard_nlls.num_jacobian_evaluations();
//     std::cout 
//        << "Testing fdf_nlls_solver with Bard function:" << std::endl
//        << "   Number of residual function calls = " << n_f_calls << std::endl
//        << "   Number of jacobian function calls = " << n_j_calls << std::endl
//        << "   Final solver status = " << gsl_strerror(slvr_status) << std::endl;

    BOOST_CHECK_EQUAL((int)slvr_status, (int)GSL_SUCCESS);
    BOOST_CHECK(n_f_calls < 10);
    BOOST_CHECK(n_j_calls <= n_f_calls);
    BOOST_CHECK( (abs(gsl_blas_dnrm2(r)-0.090636) < 0.00001)
                 || (abs(gsl_blas_dnrm2(r)-4.17476) < 0.0001) );
    if(abs(gsl_blas_dnrm2(r)-0.090636) < 0.00001) {
       BOOST_CHECK_CLOSE(gsl_vector_get(x, 0), 0.082411, 0.01);
       BOOST_CHECK_CLOSE(gsl_vector_get(x, 1), 1.1330, 0.01);
       BOOST_CHECK_CLOSE(gsl_vector_get(x, 2), 2.3437, 0.01);
    }

    if(J) {gsl_matrix_free(J); J = (gsl_matrix *)0;}
    if(r) {gsl_vector_free(r); r = (gsl_vector *)0;}
    if(x) {gsl_vector_free(x); x = (gsl_vector *)0;}
    gsl_vector_free(x0);
}



BOOST_AUTO_TEST_CASE(brown)
{
    BrownNLLSProblem brown_nlls;

    /* selection of the problem */
    f = gsl_get_lsp_fdf(&brown_nlls);

    /* for the initial guess */
    x0 = gsl_vector_calloc(f.p);

    /* comment out any of the following 3 lines if the optional 
     * solution, residual and Jacobian at the final point are not
     * desired. */
    x = gsl_vector_alloc(f.p);
    r = gsl_vector_alloc(f.n);
    /*J = gsl_matrix_alloc(f.n,f.p);*/

    /* selection of the initial guess */
    gsl_vector_set(x0, 0, 1.0);
    gsl_vector_set(x0, 1, 1.0);

    /* solve */
    slvr_status = fdf_nlls_solver(&f, x0, T, dx_epsabs, dx_epsrel, g_epsabs, 100, verbose, &num_iter, x, r, J);
    int n_f_calls = brown_nlls.num_residual_evaluations();
    int n_j_calls = brown_nlls.num_jacobian_evaluations();
//     std::cout 
//        << "Testing fdf_nlls_solver with Brown badly scaled function:" << std::endl
//        << "   Number of residual function calls = " << n_f_calls << std::endl
//        << "   Number of jacobian function calls = " << n_j_calls << std::endl
//        << "   Final solver status = " << gsl_strerror(slvr_status) << std::endl;

    BOOST_CHECK_EQUAL((int)slvr_status, (int)GSL_SUCCESS);
    BOOST_CHECK(n_f_calls < 20);
    BOOST_CHECK(n_j_calls <= n_f_calls);
    BOOST_CHECK(abs(gsl_blas_dnrm2(r)-0.0) < 0.00001);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 0), 1.0e6, 0.01);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 1), 2.0e-6, 0.01);

    if(J) {gsl_matrix_free(J); J = (gsl_matrix *)0;}
    if(r) {gsl_vector_free(r); r = (gsl_vector *)0;}
    if(x) {gsl_vector_free(x); x = (gsl_vector *)0;}
    gsl_vector_free(x0);
}



BOOST_AUTO_TEST_CASE(freudenstein_roth)
{
    FreudensteinRothNLLSProblem freudenstein_roth_nlls;

    /* selection of the problem */
    f = gsl_get_lsp_fdf(&freudenstein_roth_nlls);

    /* for the initial guess */
    x0 = gsl_vector_calloc(f.p);

    /* comment out any of the following 3 lines if the optional 
     * solution, residual and Jacobian at the final point are not
     * desired. */
    x = gsl_vector_alloc(f.p);
    r = gsl_vector_alloc(f.n);
    /*J = gsl_matrix_alloc(f.n,f.p);*/

    /* selection of the initial guess */
    gsl_vector_set(x0, 0, 0.5);
    gsl_vector_set(x0, 1, -2.0);

    /* solve */
    slvr_status = fdf_nlls_solver(&f, x0, T, dx_epsabs, dx_epsrel, g_epsabs, 100, verbose, &num_iter, x, r, J);
    int n_f_calls = freudenstein_roth_nlls.num_residual_evaluations();
    int n_j_calls = freudenstein_roth_nlls.num_jacobian_evaluations();
//     std::cout 
//        << "Testing fdf_nlls_solver with Freudenstein & Roth function:" << std::endl
//        << "   Number of residual function calls = " << n_f_calls << std::endl
//        << "   Number of jacobian function calls = " << n_j_calls << std::endl
//        << "   Final solver status = " << gsl_strerror(slvr_status) << std::endl;

    BOOST_CHECK_EQUAL((int)slvr_status, (int)GSL_SUCCESS);
    BOOST_CHECK(n_f_calls < 20);
    BOOST_CHECK(n_j_calls <= n_f_calls);
    BOOST_CHECK( (abs(gsl_blas_dnrm2(r)-0.0) < 0.00001)
                 || (abs(gsl_blas_dnrm2(r)-6.99888) < 0.0001) );
    if(abs(gsl_blas_dnrm2(r)-0.0) < 0.00001) {
       BOOST_CHECK_CLOSE(gsl_vector_get(x, 0), 5.0, 0.01);
       BOOST_CHECK_CLOSE(gsl_vector_get(x, 1), 4.0, 0.01);
    } else {
       BOOST_CHECK_CLOSE(gsl_vector_get(x, 0), 11.413, 0.01);
       BOOST_CHECK_CLOSE(gsl_vector_get(x, 1), -0.89681, 0.01);
    }

    if(J) {gsl_matrix_free(J); J = (gsl_matrix *)0;}
    if(r) {gsl_vector_free(r); r = (gsl_vector *)0;}
    if(x) {gsl_vector_free(x); x = (gsl_vector *)0;}
    gsl_vector_free(x0);
}



BOOST_AUTO_TEST_CASE(helical_valley)
{
    HelicalValleyNLLSProblem helical_valley_nlls;

    /* selection of the problem */
    f = gsl_get_lsp_fdf(&helical_valley_nlls);

    /* for the initial guess */
    x0 = gsl_vector_calloc(f.p);

    /* comment out any of the following 3 lines if the optional 
     * solution, residual and Jacobian at the final point are not
     * desired. */
    x = gsl_vector_alloc(f.p);
    r = gsl_vector_alloc(f.n);
    /*J = gsl_matrix_alloc(f.n,f.p);*/

    /* selection of the initial guess */
    gsl_vector_set(x0, 0, -1.0);
    gsl_vector_set(x0, 1, 0.0);
    gsl_vector_set(x0, 2, 0.0);

    /* solve */
    slvr_status = fdf_nlls_solver(&f, x0, T, dx_epsabs, dx_epsrel, g_epsabs, 100, verbose, &num_iter, x, r, J);
    int n_f_calls = helical_valley_nlls.num_residual_evaluations();
    int n_j_calls = helical_valley_nlls.num_jacobian_evaluations();
//     std::cout 
//        << "Testing fdf_nlls_solver with the helical valley function:" << std::endl
//        << "   Number of residual function calls = " << n_f_calls << std::endl
//        << "   Number of jacobian function calls = " << n_j_calls << std::endl
//        << "   Final solver status = " << gsl_strerror(slvr_status) << std::endl;

    BOOST_CHECK_EQUAL((int)slvr_status, (int)GSL_SUCCESS);
    BOOST_CHECK(n_f_calls < 15);
    BOOST_CHECK(n_j_calls <= n_f_calls);
    BOOST_CHECK(abs(gsl_blas_dnrm2(r)-0.0) < 0.00001);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 0), 1.0, 0.01);
    BOOST_CHECK(abs(gsl_vector_get(x, 1)-0.0) < 0.00001);
    BOOST_CHECK(abs(gsl_vector_get(x, 2)-0.0) < 0.00001);

    if(J) {gsl_matrix_free(J); J = (gsl_matrix *)0;}
    if(r) {gsl_vector_free(r); r = (gsl_vector *)0;}
    if(x) {gsl_vector_free(x); x = (gsl_vector *)0;}
    gsl_vector_free(x0);
}



BOOST_AUTO_TEST_CASE(jennrich_sampson)
{
    JennrichSampsonNLLSProblem jennrich_sampson_nlls;

    /* selection of the problem */
    f = gsl_get_lsp_fdf(&jennrich_sampson_nlls);

    /* for the initial guess */
    x0 = gsl_vector_calloc(f.p);

    /* comment out any of the following 3 lines if the optional 
     * solution, residual and Jacobian at the final point are not
     * desired. */
    x = gsl_vector_alloc(f.p);
    r = gsl_vector_alloc(f.n);
    /*J = gsl_matrix_alloc(f.n,f.p);*/

    /* selection of the initial guess */
    gsl_vector_set(x0, 0, 0.3);
    gsl_vector_set(x0, 1, 0.4);

    /* solve */
    slvr_status = fdf_nlls_solver(&f, x0, T, dx_epsabs, dx_epsrel, g_epsabs, 100, verbose, &num_iter, x, r, J);
    int n_f_calls = jennrich_sampson_nlls.num_residual_evaluations();
    int n_j_calls = jennrich_sampson_nlls.num_jacobian_evaluations();
//     std::cout 
//        << "Testing fdf_nlls_solver with Jennrich & Sampson function:" << std::endl
//        << "   Number of residual function calls = " << n_f_calls << std::endl
//        << "   Number of jacobian function calls = " << n_j_calls << std::endl
//        << "   Final solver status = " << gsl_strerror(slvr_status) << std::endl;

    BOOST_CHECK_EQUAL((int)slvr_status, (int)GSL_SUCCESS);
    BOOST_CHECK(n_f_calls < 25);
    BOOST_CHECK(n_j_calls <= n_f_calls);
    BOOST_CHECK_CLOSE(gsl_blas_dnrm2(r), 11.1518, 0.01);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 0), 0.2578, 0.01);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 1), 0.2578, 0.01);

    if(J) {gsl_matrix_free(J); J = (gsl_matrix *)0;}
    if(r) {gsl_vector_free(r); r = (gsl_vector *)0;}
    if(x) {gsl_vector_free(x); x = (gsl_vector *)0;}
    gsl_vector_free(x0);
}



BOOST_AUTO_TEST_CASE(meyer)
{
    MeyerNLLSProblem meyer_nlls;

    /* selection of the problem */
    f = gsl_get_lsp_fdf(&meyer_nlls);

    /* for the initial guess */
    x0 = gsl_vector_calloc(f.p);

    /* comment out any of the following 3 lines if the optional 
     * solution, residual and Jacobian at the final point are not
     * desired. */
    x = gsl_vector_alloc(f.p);
    r = gsl_vector_alloc(f.n);
    /*J = gsl_matrix_alloc(f.n,f.p);*/

    /* selection of the initial guess */
    gsl_vector_set(x0, 0, 0.02);
    gsl_vector_set(x0, 1, 4000);
    gsl_vector_set(x0, 2, 250);

    /* solve */
    slvr_status = fdf_nlls_solver(&f, x0, T, dx_epsabs, dx_epsrel, g_epsabs, 200, verbose, &num_iter, x, r, J);
    int n_f_calls = meyer_nlls.num_residual_evaluations();
    int n_j_calls = meyer_nlls.num_jacobian_evaluations();
//     std::cout 
//        << "Testing fdf_nlls_solver with Meyer function:" << std::endl
//        << "   Number of residual function calls = " << n_f_calls << std::endl
//        << "   Number of jacobian function calls = " << n_j_calls << std::endl
//        << "   Final solver status = " << gsl_strerror(slvr_status) << std::endl;

    BOOST_CHECK_EQUAL((int)slvr_status, (int)GSL_SUCCESS);
    BOOST_CHECK(n_f_calls < 140);
    BOOST_CHECK(n_j_calls <= n_f_calls);
    BOOST_CHECK_CLOSE(gsl_blas_dnrm2(r), 9.37794, 0.01);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 0), 0.0056096, 0.01);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 1), 6181.35, 0.01);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 2), 345.224, 0.01);

    if(J) {gsl_matrix_free(J); J = (gsl_matrix *)0;}
    if(r) {gsl_vector_free(r); r = (gsl_vector *)0;}
    if(x) {gsl_vector_free(x); x = (gsl_vector *)0;}
    gsl_vector_free(x0);
}



BOOST_AUTO_TEST_CASE(powell)
{
    PowellNLLSProblem powell_nlls;

    /* selection of the problem */
    f = gsl_get_lsp_fdf(&powell_nlls);

    /* for the initial guess */
    x0 = gsl_vector_calloc(f.p);

    /* comment out any of the following 3 lines if the optional 
     * solution, residual and Jacobian at the final point are not
     * desired. */
    x = gsl_vector_alloc(f.p);
    r = gsl_vector_alloc(f.n);
    /*J = gsl_matrix_alloc(f.n,f.p);*/

    /* selection of the initial guess */
    gsl_vector_set(x0, 0, 0.0);
    gsl_vector_set(x0, 1, 1.0);

    /* solve */
    slvr_status = fdf_nlls_solver(&f, x0, T, dx_epsabs, dx_epsrel, g_epsabs, 100, verbose, &num_iter, x, r, J);
    int n_f_calls = powell_nlls.num_residual_evaluations();
    int n_j_calls = powell_nlls.num_jacobian_evaluations();
//     std::cout 
//        << "Testing fdf_nlls_solver with Powell badly scaled function:" << std::endl
//        << "   Number of residual function calls = " << n_f_calls << std::endl
//        << "   Number of jacobian function calls = " << n_j_calls << std::endl
//        << "   Final solver status = " << gsl_strerror(slvr_status) << std::endl;

    BOOST_CHECK_EQUAL((int)slvr_status, (int)GSL_SUCCESS);
    BOOST_CHECK(n_f_calls < 20);
    BOOST_CHECK(n_j_calls <= n_f_calls);
    BOOST_CHECK(abs(gsl_blas_dnrm2(r)-0.0) < 0.00001);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 0), 1.09816e-5, 0.01);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 1), 9.10615, 0.01);

    if(J) {gsl_matrix_free(J); J = (gsl_matrix *)0;}
    if(r) {gsl_vector_free(r); r = (gsl_vector *)0;}
    if(x) {gsl_vector_free(x); x = (gsl_vector *)0;}
    gsl_vector_free(x0);
}



BOOST_AUTO_TEST_CASE(powell_singular)
{
    PowellSingularNLLSProblem powell_singular_nlls;

    /* selection of the problem */
    f = gsl_get_lsp_fdf(&powell_singular_nlls);

    /* for the initial guess */
    x0 = gsl_vector_calloc(f.p);

    /* comment out any of the following 3 lines if the optional 
     * solution, residual and Jacobian at the final point are not
     * desired. */
    x = gsl_vector_alloc(f.p);
    r = gsl_vector_alloc(f.n);
    /*J = gsl_matrix_alloc(f.n,f.p);*/

    /* selection of the initial guess */
    gsl_vector_set(x0, 0,  3.0);
    gsl_vector_set(x0, 1, -1.0);
    gsl_vector_set(x0, 2,  0.0);
    gsl_vector_set(x0, 3,  1.0);

    /* solve */
    slvr_status = fdf_nlls_solver(&f, x0, T, 1e-13, 1e-13, 1e-13, 100, verbose, &num_iter, x, r, J);
    int n_f_calls = powell_singular_nlls.num_residual_evaluations();
    int n_j_calls = powell_singular_nlls.num_jacobian_evaluations();
//     std::cout 
//        << "Testing fdf_nlls_solver with Powell singular function:" << std::endl
//        << "   Number of residual function calls = " << n_f_calls << std::endl
//        << "   Number of jacobian function calls = " << n_j_calls << std::endl
//        << "   Final solver status = " << gsl_strerror(slvr_status) << std::endl;

    BOOST_CHECK_EQUAL((int)slvr_status, (int)GSL_SUCCESS);
    BOOST_CHECK(n_f_calls < 25);
    BOOST_CHECK(n_j_calls <= n_f_calls);
    BOOST_CHECK(abs(gsl_blas_dnrm2(r)-0.0) < 0.00001);
    BOOST_CHECK(abs(gsl_vector_get(x, 0)-0.0) < 0.00001);
    BOOST_CHECK(abs(gsl_vector_get(x, 1)-0.0) < 0.00001);
    BOOST_CHECK(abs(gsl_vector_get(x, 2)-0.0) < 0.00001);
    BOOST_CHECK(abs(gsl_vector_get(x, 3)-0.0) < 0.00001);

    if(J) {gsl_matrix_free(J); J = (gsl_matrix *)0;}
    if(r) {gsl_vector_free(r); r = (gsl_vector *)0;}
    if(x) {gsl_vector_free(x); x = (gsl_vector *)0;}
    gsl_vector_free(x0);
}



BOOST_AUTO_TEST_CASE(rosenbrock2)
{
    Rosenbrock2NLLSProblem rosenbrock2_nlls;

    /* selection of the problem */
    f = gsl_get_lsp_fdf(&rosenbrock2_nlls);

    /* for the initial guess */
    x0 = gsl_vector_calloc(f.p);

    /* comment out any of the following 3 lines if the optional 
     * solution, residual and Jacobian at the final point are not
     * desired. */
    x = gsl_vector_alloc(f.p);
    r = gsl_vector_alloc(f.n);
    /*J = gsl_matrix_alloc(f.n,f.p);*/

    /* selection of the initial guess */
    gsl_vector_set(x0, 0, -1.2);
    gsl_vector_set(x0, 1, 1.0);

    /* solve */
    slvr_status = fdf_nlls_solver(&f, x0, T, dx_epsabs, dx_epsrel, g_epsabs, 100, verbose, &num_iter, x, r, J);
    int n_f_calls = rosenbrock2_nlls.num_residual_evaluations();
    int n_j_calls = rosenbrock2_nlls.num_jacobian_evaluations();
//     std::cout 
//        << "Testing fdf_nlls_solver with Rosenbrock function:" << std::endl
//        << "   Number of residual function calls = " << n_f_calls << std::endl
//        << "   Number of jacobian function calls = " << n_j_calls << std::endl
//        << "   Final solver status = " << gsl_strerror(slvr_status) << std::endl;

    BOOST_CHECK_EQUAL((int)slvr_status, (int)GSL_SUCCESS);
    BOOST_CHECK(n_f_calls < 25);
    BOOST_CHECK(n_j_calls <= n_f_calls);
    BOOST_CHECK(abs(gsl_blas_dnrm2(r)-0.0) < 0.00001);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 0), 1.0, 0.01);
    BOOST_CHECK_CLOSE(gsl_vector_get(x, 1), 1.0, 0.01);

    if(J) {gsl_matrix_free(J); J = (gsl_matrix *)0;}
    if(r) {gsl_vector_free(r); r = (gsl_vector *)0;}
    if(x) {gsl_vector_free(x); x = (gsl_vector *)0;}
    gsl_vector_free(x0);
}



BOOST_AUTO_TEST_SUITE_END()
