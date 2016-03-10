#include <unit_test_support.h>
#include <fp_exception.h>
#include <nlls_max_a_posteriori.h>

using namespace FullPhysics;
using namespace blitz;

class MaxAPosterioriTest : public MaxAPosteriori
{
public:
  MaxAPosterioriTest(const Array<double, 1> msrmnt,
                     const Array<double, 1> msrmnt_err_cov,
                     const Array<double, 1> a_priori_params,
                     const Array<double, 2> a_priori_cov,
                     const Array<double, 1> mdl,
                     const Array<double, 2> mdl_jac)
    : ModelMeasure(msrmnt, msrmnt_err_cov), MaxAPosteriori(a_priori_params, a_priori_cov),
      Mdl(mdl.copy()), Mdl_jac(mdl_jac.copy())
  {}

  virtual ~MaxAPosterioriTest() {}
  virtual void model_eval() {M.reference(Mdl);}
  virtual void jacobian_eval() {K.reference(Mdl_jac);}
  virtual int expected_parameter_size() const { return 3; }
private:
  Array<double, 1> Mdl;
  Array<double, 2> Mdl_jac;
};


BOOST_FIXTURE_TEST_SUITE(nlls_max_a_posteriori, GlobalFixture)


BOOST_AUTO_TEST_CASE(residual_jacobian)
{
  Array<double, 1> Msrmnt(4);
  Array<double, 1> Msrmnt_err_cov(4);
  Array<double, 1> A_p_params(3);
  Array<double, 2> A_p_cov(3,3);
  Array<double, 2> A_p_cov_chol_inv(3,3);
  Array<double, 1> Mdl(4);
  Array<double, 2> Mdl_jac(4,3);
  Array<double, 1> another_mdl;
  Array<double, 2> another_mdl_jac;

  Array<double, 1> Params(3);

  Msrmnt = 1.5, 2.5, 0.75, -0.5;
  Msrmnt_err_cov = 4.0, 9.0, 1.0, 16.0;

  A_p_params = -0.5, 0.75, 0.095;

  A_p_cov =
    2.2500,     3.00000,    0.750000,
    3.0000,   114.25000,   40.375000,
    0.7500,    40.37500,   39.312500;

  // Inverse of the Cholesky decomposition of the above matrix
  A_p_cov_chol_inv =
    0.666666666666667,  0.000000000000000, -0.000000000000000,
   -0.126984126984127,  0.095238095238095, -0.000000000000000,
    0.028571428571429, -0.071428571428571,  0.200000000000000;


  Mdl = 1.25, 2.75, 0.25, -0.25;

  Mdl_jac = 
     1.00,  2.00,  3.50,
    -0.50,  0.25, 10.50,
     2.50, -0.75,  1.25,
     3.00, -3.15, -4.35;

  Params = 0.5, 1.75, 1.50;


  Array<double, 1> mle_res(4);
  Array<double, 1> ap_res(3);
  Array<double, 1> res(7);
  Array<double, 2> mle_jac(4,3);
  Array<double, 2> ap_jac(3,3);
  Array<double, 2> jac(7,3);
  Array<double, 1> grad(3);
  Array<double, 1> another_res;
  Array<double, 2> another_jac;
  Array<double, 2> a_post_cov(3,3);
  Array<double, 2> avg_krnl(3,3);
  Array<double, 1> x_a_prior_uncert(3);

  mle_res = (1.25-1.5)/2.0,   (2.75-2.5)/3.0,   (0.25-0.75)/1.0,   (-0.25+0.5)/4.0;

  ap_res = 
    (0.5-(-0.5))*0.666666666666667 + (1.75-0.75)*0.000000000000000 + (1.50-0.095)*0.000000000000000, 
    (0.5-(-0.5))*(-0.126984126984127) + (1.75-0.75)*0.095238095238095 + (1.50-0.095)*0.000000000000000,
    (0.5-(-0.5))*0.028571428571429 + (1.75-0.75)*(-0.071428571428571) + (1.50-0.095)*0.200000000000000;

  res = // mle_res augmented by ap_res
    (1.25-1.5)/2.0,   (2.75-2.5)/3.0,   (0.25-0.75)/1.0,   (-0.25+0.5)/4.0,
    (0.5-(-0.5))*0.666666666666667 + (1.75-0.75)*0.000000000000000 + (1.50-0.095)*0.000000000000000, 
    (0.5-(-0.5))*(-0.126984126984127) + (1.75-0.75)*0.095238095238095 + (1.50-0.095)*0.000000000000000,
    (0.5-(-0.5))*0.028571428571429 + (1.75-0.75)*(-0.071428571428571) + (1.50-0.095)*0.200000000000000;

  mle_jac = 
     1.00/2.0,  2.00/2.0,  3.50/2.0,
    -0.50/3.0,  0.25/3.0, 10.50/3.0,
     2.50/1.0, -0.75/1.0,  1.25/1.0,
     3.00/4.0, -3.15/4.0, -4.35/4.0;

  // a-posteriori jacobian portion must be A_p_cov_chol_inv

  jac = // mle_jac augmented by A_p_cov_chol_inv (the a-posteriori portion of the jacobian)
     1.00/2.0,  2.00/2.0,  3.50/2.0,
    -0.50/3.0,  0.25/3.0, 10.50/3.0,
     2.50/1.0, -0.75/1.0,  1.25/1.0,
     3.00/4.0, -3.15/4.0, -4.35/4.0,
     0.666666666666667,  0.000000000000000, -0.000000000000000,
    -0.126984126984127,  0.095238095238095, -0.000000000000000,
     0.028571428571429, -0.071428571428571,  0.200000000000000;


  grad =
    ((1.25-1.5)/2.0)*(1.00/2.0) + ((2.75-2.5)/3.0)*(-0.50/3.0) + ((0.25-0.75)/1.0)*(2.50/1.0) + ((-0.25+0.5)/4.0)*(3.00/4.0)
       + res(4)*0.666666666666667 + res(5)*(-0.126984126984127) + res(6)*0.028571428571429,
    ((1.25-1.5)/2.0)*(2.00/2.0) + ((2.75-2.5)/3.0)*(0.25/3.0) + ((0.25-0.75)/1.0)*(-0.75/1.0) + ((-0.25+0.5)/4.0)*(-3.15/4.0)
       + res(4)*0.000000000000000 + res(5)*0.095238095238095 + res(6)*(-0.071428571428571),
    ((1.25-1.5)/2.0)*(3.50/2.0) + ((2.75-2.5)/3.0)*(10.50/3.0) + ((0.25-0.75)/1.0)*(1.25/1.0) + ((-0.25+0.5)/4.0)*(-4.35/4.0)
       + res(4)*0.000000000000000 + res(5)*0.000000000000000 + res(6)*0.200000000000000;

  a_post_cov =  // hand calculated
    0.218174,  0.248751, -0.058177, 
    0.248751,  0.785002, -0.120252, 
   -0.058177, -0.120252, 0.0765678;

  avg_krnl = // hand calculated
    0.90318596,   -0.0012727012,   0.0046342441,
   -0.10298767,    0.99067278,     0.014603333,
    0.024704907,   0.0019757757,   0.99555177;

  x_a_prior_uncert = sqrt(2.2500), sqrt(114.2500), sqrt(39.312500);

  boost::shared_ptr<MaxAPosterioriTest> mapt(new MaxAPosterioriTest(Msrmnt, Msrmnt_err_cov, A_p_params, A_p_cov, Mdl, Mdl_jac));
  NLLSMaxAPosteriori nlls_map(mapt);
  nlls_map.parameters(Params);

  mapt->model_jacobian_x(Params, another_mdl, another_mdl_jac);

  BOOST_CHECK(sum(abs(mapt->measurement()-Msrmnt)) < 0.000001);
  BOOST_CHECK(sum(abs(mapt->measurement_error_cov()-Msrmnt_err_cov)) < 0.000001);
  BOOST_CHECK(sum(abs(mapt->a_priori_params()-A_p_params)) < 0.000001);
  BOOST_CHECK(sum(abs(mapt->a_priori_cov()-A_p_cov)) < 0.000001);
  BOOST_CHECK(sum(abs(mapt->model()-Mdl)) < 0.000001);
  BOOST_CHECK(sum(abs(mapt->jacobian()-Mdl_jac)) < 0.000001);
  BOOST_CHECK_EQUAL(mapt->measurement_size(), 4);
  BOOST_CHECK_EQUAL(mapt->parameter_size(), 3);
  BOOST_CHECK(sum(abs(mapt->a_posteriori_covariance()-a_post_cov)) < 0.00001);
  BOOST_CHECK(sum(abs(mapt->averaging_kernel()-avg_krnl)) < 0.00001);
  BOOST_CHECK(sum(abs(mapt->param_a_priori_uncertainty()-x_a_prior_uncert)) < 0.000001);

  BOOST_CHECK_EQUAL(nlls_map.residual_size(), 7);
  BOOST_CHECK_EQUAL(nlls_map.expected_parameter_size(), 3);
  BOOST_CHECK_EQUAL(nlls_map.parameter_size(), 3);
  //
  BOOST_CHECK_EQUAL(nlls_map.num_residual_evaluations(), 0);
  BOOST_CHECK_EQUAL(nlls_map.num_jacobian_evaluations(), 0);
  BOOST_CHECK(sum(abs(nlls_map.residual()-res)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_map.num_jacobian_evaluations(), 0);
  BOOST_CHECK(sum(abs(nlls_map.jacobian()-jac)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_map.num_jacobian_evaluations(), 1);
  BOOST_CHECK(sum(abs(nlls_map.gradient()-grad)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_map.num_jacobian_evaluations(), 1);
  nlls_map.residual_jacobian_x(Params, another_res, another_jac);
  BOOST_CHECK(sum(abs(another_jac-jac)) < 0.000001);
  BOOST_CHECK(sum(abs(another_res-res)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_map.num_jacobian_evaluations(), 1);
  Params(0) = -10.0;  // This changes the residual and gradient
  nlls_map.residual_jacobian_x(Params, another_res, another_jac);
  BOOST_CHECK(sum(abs(another_jac-jac)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map.num_residual_evaluations(), 2);
  BOOST_CHECK_EQUAL(nlls_map.num_jacobian_evaluations(), 2);
  Params(0) = 3.0;
  nlls_map.parameters(Params);
  BOOST_CHECK(sum(abs(nlls_map.jacobian()-jac)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map.num_residual_evaluations(), 2);
  BOOST_CHECK_EQUAL(nlls_map.num_jacobian_evaluations(), 3);

  NLLSMaxAPosteriori nlls_map_2(mapt, true);
  Params(0) = 0.5;  // Restore the original value of Param(0)
  nlls_map_2.parameters(Params);

  BOOST_CHECK_EQUAL(nlls_map_2.residual_size(), 7);
  BOOST_CHECK_EQUAL(nlls_map_2.expected_parameter_size(), 3);
  BOOST_CHECK_EQUAL(nlls_map_2.parameter_size(), 3);
  //
  BOOST_CHECK_EQUAL(nlls_map_2.num_residual_evaluations(), 0);
  BOOST_CHECK_EQUAL(nlls_map_2.num_jacobian_evaluations(), 0);
  BOOST_CHECK(sum(abs(nlls_map_2.residual()-res)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map_2.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_map_2.num_jacobian_evaluations(), 1);
  BOOST_CHECK(sum(abs(nlls_map_2.jacobian()-jac)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map_2.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_map_2.num_jacobian_evaluations(), 1);
  BOOST_CHECK(sum(abs(nlls_map_2.gradient()-grad)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map_2.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_map_2.num_jacobian_evaluations(), 1);
  nlls_map_2.residual_jacobian_x(Params, another_res, another_jac);
  BOOST_CHECK(sum(abs(another_jac-jac)) < 0.000001);
  BOOST_CHECK(sum(abs(another_res-res)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map_2.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_map_2.num_jacobian_evaluations(), 1);
  Params(0) = 5.0;  // This changes the residual and gradient
  nlls_map_2.residual_jacobian_x(Params, another_res, another_jac);
  BOOST_CHECK(sum(abs(another_jac-jac)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map_2.num_residual_evaluations(), 2);
  BOOST_CHECK_EQUAL(nlls_map_2.num_jacobian_evaluations(), 2);
  Params(0) = 3.0;  // This changes the residual and gradient
  nlls_map_2.parameters(Params);
  BOOST_CHECK(sum(abs(nlls_map_2.jacobian()-jac)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_map_2.num_residual_evaluations(), 3);
  BOOST_CHECK_EQUAL(nlls_map_2.num_jacobian_evaluations(), 3);

}


BOOST_AUTO_TEST_SUITE_END()
