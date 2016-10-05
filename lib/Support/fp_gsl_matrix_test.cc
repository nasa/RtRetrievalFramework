#include "fp_gsl_matrix.h"
#include "linear_algebra.h"
#include "unit_test_support.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(fp_gsl_matrix, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  blitz::Array<double, 2> t1(2, 3);
  t1 = 
    1, 2, 3,
    5, 5, 6;
  GslMatrix m1(t1);
  BOOST_CHECK(blitz::all(fabs(m1.blitz_array() - t1) < 1e-8));
  BOOST_CHECK_EQUAL(m1.gsl()->size1, (size_t) 2);
  BOOST_CHECK_EQUAL(m1.gsl()->size2, (size_t) 3);
  for(int i = 0; i < t1.rows(); ++i)
    for(int j = 0; j < t1.cols(); ++j)
      BOOST_CHECK_CLOSE(gsl_matrix_get(m1.gsl(), i, j), t1(i, j), 1e-4);
  gsl_matrix_set(m1.gsl(), 1, 2, 3.0);
  BOOST_CHECK_CLOSE(t1(1, 2), 3.0, 1e-4);
  gsl_matrix* t2 = gsl_matrix_alloc(2, 3);
  GslMatrix m2(t2);
  m2.blitz_array() = t1;
  BOOST_CHECK(blitz::all(fabs(m2.blitz_array() - t1) < 1e-8));
  BOOST_CHECK_EQUAL(m2.gsl()->size1, (size_t) 2);
  BOOST_CHECK_EQUAL(m2.gsl()->size2, (size_t) 3);
  for(int i = 0; i < t1.rows(); ++i)
    for(int j = 0; j < t1.cols(); ++j)
      BOOST_CHECK_CLOSE(gsl_matrix_get(m2.gsl(), i, j), t1(i, j), 1e-4);
}

BOOST_AUTO_TEST_CASE(fortran_order)
{
  blitz::Array<double, 2> t1_c(2, 3);
  t1_c = 
    1, 2, 3,
    5, 5, 6;
  blitz::Array<double, 2> t1(to_fortran(t1_c));
  GslMatrix m1(t1);
  BOOST_CHECK(blitz::all(fabs(m1.blitz_array() - t1) < 1e-8));
  BOOST_CHECK_EQUAL(m1.gsl()->size1, (size_t) 2);
  BOOST_CHECK_EQUAL(m1.gsl()->size2, (size_t) 3);
  for(int i = 0; i < t1.rows(); ++i)
    for(int j = 0; j < t1.cols(); ++j)
      BOOST_CHECK_CLOSE(gsl_matrix_get(m1.gsl(), i, j), t1(i, j), 1e-4);
  gsl_matrix_set(m1.gsl(), 1, 2, 3.0);
  gsl_matrix* t2 = gsl_matrix_alloc(2, 3);
  GslMatrix m2(t2);
  m2.blitz_array() = t1;
  BOOST_CHECK(blitz::all(fabs(m2.blitz_array() - t1) < 1e-8));
  BOOST_CHECK_EQUAL(m2.gsl()->size1, (size_t) 2);
  BOOST_CHECK_EQUAL(m2.gsl()->size2, (size_t) 3);
  for(int i = 0; i < t1.rows(); ++i)
    for(int j = 0; j < t1.cols(); ++j)
      BOOST_CHECK_CLOSE(gsl_matrix_get(m2.gsl(), i, j), t1(i, j), 1e-4);
}

BOOST_AUTO_TEST_CASE(fp_gsl_vector)
{
  blitz::Array<double, 1> t1(3);
  t1 = 1, 2, 3;
  GslVector m1(t1);
  BOOST_CHECK(blitz::all(fabs(m1.blitz_array() - t1) < 1e-8));
  BOOST_CHECK_EQUAL(m1.gsl()->size, (size_t) 3);
  for(int i = 0; i < t1.rows(); ++i)
    BOOST_CHECK_CLOSE(gsl_vector_get(m1.gsl(), i), t1(i), 1e-4);
  gsl_vector_set(m1.gsl(), 1, 3.0);
  BOOST_CHECK_CLOSE(t1(1), 3.0, 1e-4);
  gsl_vector* t2 = gsl_vector_alloc(3);
  GslVector m2(t2);
  m2.blitz_array() = t1;
  BOOST_CHECK(blitz::all(fabs(m2.blitz_array() - t1) < 1e-8));
  BOOST_CHECK_EQUAL(m2.gsl()->size, (size_t) 3);
  for(int i = 0; i < t1.rows(); ++i)
    BOOST_CHECK_CLOSE(gsl_vector_get(m2.gsl(), i), t1(i), 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
