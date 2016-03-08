#include <unit_test_support.h>
#include <fp_exception.h>
#include <nlls_max_likelihood.h>

using namespace FullPhysics;
using namespace blitz;

class MaxLikelihoodTest : public MaxLikelihood
{
public:
  MaxLikelihoodTest(const Array<double, 1> msrmnt, 
                    const Array<double, 1> msrmnt_err_cov, 
                    const Array<double, 1> mdl, 
                    const Array<double, 2> mdl_jac)
    : ModelMeasure(msrmnt, msrmnt_err_cov), Mdl(mdl.copy()), Mdl_jac(mdl_jac.copy())
  {}

  virtual ~MaxLikelihoodTest() {}
  virtual void model_eval() {M.reference(Mdl);}
  virtual void jacobian_eval() {K.reference(Mdl_jac);}
  virtual int expected_parameter_size() const { return 3; }
private:
  Array<double, 1> Mdl;
  Array<double, 2> Mdl_jac;
};


BOOST_FIXTURE_TEST_SUITE(nlls_max_likelihood, GlobalFixture)


BOOST_AUTO_TEST_CASE(residual_jacobian)
{
  Array<double, 1> Msrmnt(4);
  Array<double, 1> Msrmnt_err_cov(4);
  Array<double, 1> Mdl(4);
  Array<double, 2> Mdl_jac(4,3);
  Array<double, 1> another_mdl;
  Array<double, 2> another_mdl_jac;

  Array<double, 1> Params(3);

  Msrmnt = 1.5, 2.5, 0.75, -0.5;
  Msrmnt_err_cov = 4.0, 9.0, 1.0, 16.0;
  Mdl = 1.25, 2.75, 0.25, -0.25;

  Mdl_jac = 
     1.00,  2.00,  3.50,
    -0.50,  0.25, 10.50,
     2.50, -0.75,  1.25,
     3.00, -3.15, -4.35;

  Params = 0.5, 1.75, 1.50;

  Array<double, 1> res(4);
  Array<double, 2> jac(4,3);
  Array<double, 1> grad(3);
  Array<double, 1> another_res;
  Array<double, 2> another_jac;

  res = (1.25-1.5)/2.0, (2.75-2.5)/3.0, (0.25-0.75)/1.0, (-0.25+0.5)/4.0;

  jac = 
     1.00/2.0,  2.00/2.0,  3.50/2.0,
    -0.50/3.0,  0.25/3.0, 10.50/3.0,
     2.50/1.0, -0.75/1.0,  1.25/1.0,
     3.00/4.0, -3.15/4.0, -4.35/4.0;

  grad =
    ((1.25-1.5)/2.0)*(1.00/2.0) + ((2.75-2.5)/3.0)*(-0.50/3.0) + ((0.25-0.75)/1.0)*(2.50/1.0) + ((-0.25+0.5)/4.0)*(3.00/4.0),
    ((1.25-1.5)/2.0)*(2.00/2.0) + ((2.75-2.5)/3.0)*(0.25/3.0) + ((0.25-0.75)/1.0)*(-0.75/1.0) + ((-0.25+0.5)/4.0)*(-3.15/4.0),
    ((1.25-1.5)/2.0)*(3.50/2.0) + ((2.75-2.5)/3.0)*(10.50/3.0) + ((0.25-0.75)/1.0)*(1.25/1.0) + ((-0.25+0.5)/4.0)*(-4.35/4.0);

  boost::shared_ptr<MaxLikelihoodTest> mlt(new MaxLikelihoodTest(Msrmnt, Msrmnt_err_cov, Mdl, Mdl_jac));
  NLLSMaxLikelihood nlls_ml(mlt);
  nlls_ml.parameters(Params);

  mlt->model_jacobian_x(Params, another_mdl, another_mdl_jac);

  BOOST_CHECK(sum(abs(mlt->measurement()-Msrmnt)) < 0.000001);
  BOOST_CHECK(sum(abs(mlt->measurement_error_cov()-Msrmnt_err_cov)) < 0.000001);
  BOOST_CHECK(sum(abs(mlt->model()-Mdl)) < 0.000001);
  BOOST_CHECK(sum(abs(mlt->jacobian()-Mdl_jac)) < 0.000001);
  BOOST_CHECK(sum(abs(another_mdl-Mdl)) < 0.000001);
  BOOST_CHECK(sum(abs(another_mdl_jac-Mdl_jac)) < 0.000001);
  BOOST_CHECK_EQUAL(mlt->measurement_size(), 4);
  BOOST_CHECK_EQUAL(mlt->parameter_size(), 3);
  BOOST_CHECK_MATRIX_CLOSE(mlt->parameters(), Params);

  BOOST_CHECK_EQUAL(nlls_ml.residual_size(), 4);
  BOOST_CHECK_EQUAL(nlls_ml.expected_parameter_size(), 3);
  BOOST_CHECK_EQUAL(nlls_ml.parameter_size(), 3);
  //
  BOOST_CHECK_EQUAL(nlls_ml.num_residual_evaluations(), 0);
  BOOST_CHECK_EQUAL(nlls_ml.num_jacobian_evaluations(), 0);
  BOOST_CHECK(sum(abs(nlls_ml.residual()-res)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_ml.num_jacobian_evaluations(), 0);
  BOOST_CHECK(sum(abs(nlls_ml.jacobian()-jac)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_ml.num_jacobian_evaluations(), 1);
  BOOST_CHECK(sum(abs(nlls_ml.gradient()-grad)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_ml.num_jacobian_evaluations(), 1);
  nlls_ml.residual_jacobian_x(Params, another_res, another_jac);
  BOOST_CHECK(sum(abs(another_jac-jac)) < 0.000001);
  BOOST_CHECK(sum(abs(another_res-res)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_ml.num_jacobian_evaluations(), 1);
  Params(0) = -10.0;
  nlls_ml.residual_jacobian_x(Params, another_res, another_jac);
  BOOST_CHECK(sum(abs(another_jac-jac)) < 0.000001);
  BOOST_CHECK(sum(abs(another_res-res)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml.num_residual_evaluations(), 2);
  BOOST_CHECK_EQUAL(nlls_ml.num_jacobian_evaluations(), 2);
  Params(0) = 3.0;
  nlls_ml.parameters(Params);
  BOOST_CHECK(sum(abs(nlls_ml.jacobian()-jac)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml.num_residual_evaluations(), 2);
  BOOST_CHECK_EQUAL(nlls_ml.num_jacobian_evaluations(), 3);

  NLLSMaxLikelihood nlls_ml_2(mlt, true);
  nlls_ml_2.parameters(Params);

  BOOST_CHECK_EQUAL(nlls_ml_2.residual_size(), 4);
  BOOST_CHECK_EQUAL(nlls_ml_2.expected_parameter_size(), 3);
  BOOST_CHECK_EQUAL(nlls_ml_2.parameter_size(), 3);
  //
  BOOST_CHECK_EQUAL(nlls_ml_2.num_residual_evaluations(), 0);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_jacobian_evaluations(), 0);
  BOOST_CHECK(sum(abs(nlls_ml_2.residual()-res)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_jacobian_evaluations(), 1);
  BOOST_CHECK(sum(abs(nlls_ml_2.jacobian()-jac)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_jacobian_evaluations(), 1);
  BOOST_CHECK(sum(abs(nlls_ml_2.gradient()-grad)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_jacobian_evaluations(), 1);
  nlls_ml_2.residual_jacobian_x(Params, another_res, another_jac);
  BOOST_CHECK(sum(abs(another_jac-jac)) < 0.000001);
  BOOST_CHECK(sum(abs(another_res-res)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_residual_evaluations(), 1);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_jacobian_evaluations(), 1);
  Params(0) = 5.0;
  nlls_ml_2.residual_jacobian_x(Params, another_res, another_jac);
  BOOST_CHECK(sum(abs(another_jac-jac)) < 0.000001);
  BOOST_CHECK(sum(abs(another_res-res)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_residual_evaluations(), 2);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_jacobian_evaluations(), 2);
  Params(0) = 3.0;
  nlls_ml_2.parameters(Params);
  BOOST_CHECK(sum(abs(nlls_ml_2.jacobian()-jac)) < 0.000001);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_residual_evaluations(), 3);
  BOOST_CHECK_EQUAL(nlls_ml_2.num_jacobian_evaluations(), 3);
  
}


BOOST_AUTO_TEST_SUITE_END()
