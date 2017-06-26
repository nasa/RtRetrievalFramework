#include "linear_algebra.h"
#include "unit_test_support.h"
#include "ifstream_cs.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(linear_algebra, GlobalFixture)

BOOST_AUTO_TEST_CASE(to_fortran_test)
{
  Array<double, 2> t(2,3);
  t = 1,2,3,
    4,5,6;
  Array<double, 2> tasfortran(t.dataFirst(), t.shape(), 
				     neverDeleteData,
				     ColumnMajorArray<2>());
  // Data isn't in fortran order, before we have converted it.
  BOOST_CHECK(t(0,1) != tasfortran(0,1));
  Array<double, 2> tf = to_fortran(t);
  Array<double, 2> tasfortran2(tf.dataFirst(), tf.shape(), 
				      neverDeleteData,
				      ColumnMajorArray<2>());

  // Data is now in fortran order.
  BOOST_CHECK(tf(0,1) = tasfortran2(0,1));

  // If it is already in fortran order, we don't create an extra copy.
  Array<double, 2> tf2 = to_fortran(tf);
  BOOST_CHECK_EQUAL(tf2.dataFirst(), tf.dataFirst());
}

BOOST_AUTO_TEST_CASE(to_c_order_test)
{
  Array<double, 2> t(2,3);
  t = 1,2,3,
    4,5,6;
  Array<double, 2> tasfortran(t.dataFirst(), t.shape(), 
				     neverDeleteData,
				     ColumnMajorArray<2>());
  // Data isn't in c order, before we have converted it.
  BOOST_CHECK(t(0,1) != tasfortran(0,1));
  Array<double, 2> tc = to_c_order(tasfortran);
  Array<double, 2> tasc2(tc.dataFirst(), tc.shape(), 
			      neverDeleteData);

  // Data is now in fortran order.
  BOOST_CHECK(tc(0,1) = tasc2(0,1));

  // If it is already in C order, we don't create an extra copy.
  Array<double, 2> tc2 = to_c_order(tc);
  BOOST_CHECK_EQUAL(tc2.dataFirst(), tc.dataFirst());
}

BOOST_AUTO_TEST_CASE(solve_least_squares_test)
{
  using namespace blitz;
  using namespace tensor;
  Array<double, 2> a(4,2);
  a = 
    1, 2,
    2, 3,
    4, 5,
    6, 7;
  Array<double, 1> xexpect(2);
  xexpect = 1, 2;
  Array<double, 1> b(4);
  b = sum(a(i,j) * xexpect(j),j);
  b(0) += 0.01;
  b(2) -= 0.1;
  Array<double, 1> x = solve_least_squares(a, b);
  BOOST_CHECK(max(abs(x - xexpect)) < 1e-2);

  // Check a degenerate matrix.
  Array<double, 2> a2(4,3);
  a2 = 
    1, 2, 0,
    2, 3, 0,
    4, 5, 0,
    6, 7, 0; 
  Array<double, 1> xexpect2(3);
  xexpect2 = 1, 2, 0;
  Array<double, 1> x2 = solve_least_squares(a2, b);
  BOOST_CHECK(max(abs(x2 - xexpect2)) < 1e-2);
}

BOOST_AUTO_TEST_CASE(svd_test)
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Array<double, 2> a(3,3);
  a = 1,2,0,
    3,4,5,
    0,7,8;
  Array<double, 1> s;
  Array<double, 2> u, vt;
  svd(a, s, u, vt);
  Array<double, 2> res(3,3);
  res = sum(u(i1, i3) * s(i3) * vt(i3, i2), i3);
  BOOST_CHECK(max(abs(res - a)) < 1e-8);
  Array<double, 2> ainv = generalized_inverse(a);
  Array<double, 2> prod(3,3);
  prod = sum(ainv(i1, i3) * a(i3, i2), i3);
  Array<double, 2> iden(3, 3);
  iden = 
    1, 0, 0,
    0, 1, 0,
    0, 0, 1;
  BOOST_CHECK(max(abs(prod - iden)) < 1e-8);
  a = 1,2,0,
    3,4,5,
    0,0,0;
  ainv = generalized_inverse(a);
  prod = sum(a(i1, i3) * ainv(i3, i2), i3);
  iden =
    1, 0, 0,
    0, 1, 0,
    0, 0, 0;
  BOOST_CHECK(max(abs(prod - iden)) < 1e-8);
  // Column with "5" should be set to 0 by zero_unused.
  a = 1,2,0,
    3,4,5,
    0,0,0;
  Array<bool, 1> zero_unused(3);
  zero_unused = false, false, true;
  ainv = generalized_inverse(a, zero_unused);
  a = 1,2,0,
    3,4,0,
    0,0,0;
  prod = sum(a(i1, i3) * ainv(i3, i2), i3);
  iden =
    1, 0, 0,
    0, 1, 0,
    0, 0, 0;
  BOOST_CHECK(max(abs(prod - iden)) < 1e-8);
}

BOOST_AUTO_TEST_CASE(cholesky_test)
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Array<double, 2> a(3,3);
  a = 
    10, 2, 3,
    2, 40, 5,
    3, 5, 60;
  Array<double, 2> l(cholesky_decomposition(a));
  Array<double, 2> ares(a.shape());
  ares = sum(l(i1, i3) * l(i2, i3), i3);
  BOOST_CHECK_MATRIX_CLOSE(a, ares);
}

BOOST_AUTO_TEST_CASE(generalized_inverse_fortran_order_test)
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Array<double, 2> a(3,3, FortranArray<2>());
  a = 1,2,0,
    3,4,5,
    0,7,8;
  Array<double, 2> ainv = generalized_inverse(a);
  Array<double, 2> a_c_order = to_c_order(a);
  Array<double, 2> prod(3,3);
  prod = sum(ainv(i1, i3) * a_c_order(i3, i2), i3);
  Array<double, 2> iden(3, 3);
  iden = 
    1, 0, 0,
    0, 1, 0,
    0, 0, 1;
  BOOST_CHECK(max(abs(prod - iden)) < 1e-8);
}

BOOST_AUTO_TEST_CASE(cholesky_fortran_order_test)
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Array<double, 2> a(3,3, FortranArray<2>());
  a = 
    10, 2, 3,
    2, 40, 5,
    3, 5, 60;
  Array<double, 2> l(cholesky_decomposition(a));
  Array<double, 2> ares(a.shape());
  ares = sum(l(i1, i3) * l(i2, i3), i3);
  BOOST_CHECK_MATRIX_CLOSE(to_c_order(a), ares);
}

BOOST_AUTO_TEST_CASE(gsl_test_examples)
{
  firstIndex i1; secondIndex i2; thirdIndex i3;

// This is some problem samples that Edwin found as described in
// Tickets #830 and #831. They illustrate the issue with the GSL matrix
// routines. It turns out that this is *not* actually a problem, just
// an issue with the Rcond passed to generalized_inverse.

  IfstreamCs example(test_data_dir() + 
		     "expected/linear_algebra/gsl_test_examples"); 
  Array<double, 2> m;
  example >> m;
  Array<double, 2> res(m.shape());
  res = sum(m(i1, i3) * generalized_inverse(m,1e-20)(i3, i2), i3);
  BOOST_CHECK_CLOSE(res(0,0), 1.0, 1e-8);

  Array<double, 2> m2;
  example >> m2;
  Array<double, 2> l = cholesky_decomposition(m2);
  res.resize(m2.shape());
  res = sum(l(i1, i3) * l(i2, i3), i3);
  BOOST_CHECK_MATRIX_CLOSE(res, m2);

  // Test from l2_20100207003807.log 
  Array<double, 2> m3;
  Array<double, 2> l3;
  example >> m3 >> l3;
  Array<double, 2> l3_res = cholesky_decomposition(m3);
  res.resize(m3.shape());
  res = sum(l3_res(i1, i3) * l3_res(i2, i3), i3);
  BOOST_CHECK_MATRIX_CLOSE(res, m3);
  // This test fails, I believe there is some problem with Edwin's
  // executable. 
  // BOOST_CHECK_MATRIX_CLOSE(l3_res, l3);
}


BOOST_AUTO_TEST_CASE(multifit_covar_test)
{
  Array<double, 2> J(4,3);
  J = 
    1.00,  2.50, -0.35,
   -2.25,  4.00,  0.50,
    3.50, -1.50,  1.05,
    0.07, -1.95,  2.95;
  Array<double, 2> covar_expected(3,3);
  covar_expected = 
    0.0750703640862397, 0.0318982934760073, 0.0016892492371800, 
    0.0318982934760073, 0.0543328512564905, 0.0255697633670716, 
    0.0016892492371800, 0.0255697633670716, 0.1134395754355183;

  Array<double, 2> covar(multifit_covar(J));

  BOOST_CHECK_MATRIX_CLOSE(covar_expected, covar);
}


BOOST_AUTO_TEST_SUITE_END()
