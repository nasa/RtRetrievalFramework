#ifndef UNIT_TEST_SUPPORT_H
#define UNIT_TEST_SUPPORT_H
#include "global_fixture.h"
#include "ifstream_cs.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <blitz/array.h>
#include <string>

// ***********************************************************************
// This is some common routines used in unit tests.
//
// See also global_fixture.h, which defines a global fixture that sets up
// some useful functions available everywhere.
// ***********************************************************************

//-----------------------------------------------------------------------
/// Check for two matrixes being equal
//-----------------------------------------------------------------------

#define BOOST_CHECK_MATRIX_CLOSE_TOL( A, B, tol ) BOOST_CHECK_LT(max(abs(A - (B))), tol)
#define BOOST_CHECK_MATRIX_CLOSE( A, B ) BOOST_CHECK_MATRIX_CLOSE_TOL(A, B, 1e-8)

//-----------------------------------------------------------------------
/// Check for two ArrayAd being equal
//-----------------------------------------------------------------------

#define BOOST_CHECK_ARRAYAD_CLOSE_TOL( A, B, tol ) BOOST_CHECK_MATRIX_CLOSE_TOL(A.value(), B.value(), tol); BOOST_CHECK_MATRIX_CLOSE_TOL(A.jacobian(), B.jacobian(), tol)
#define BOOST_CHECK_ARRAYAD_CLOSE( A, B ) BOOST_CHECK_MATRIX_CLOSE(A.value(), B.value()); BOOST_CHECK_MATRIX_CLOSE(A.jacobian(), B.jacobian())

//-----------------------------------------------------------------------
/// Throw Exception, with string the given value.
//-----------------------------------------------------------------------

#define CHECK_THROW_EXCEPTION( S , V ) \
  try { \
    S; \
    BOOST_ERROR("Exception was not thrown"); \
  } catch(const FullPhysics::Exception& E) {	\
    BOOST_CHECK_EQUAL(std::string(E.what()), std::string(V)) ;	\
  }

template<class T, int N>
inline void compare_shapes(const blitz::Array<T, N>& A, const blitz::Array<T, N>& B) {
    if( product(A.shape()) != product(B.shape()) ) {
        std::cerr << "Shapes of arrays do not match: " << A.shape() << " vs. " << B.shape() << std::endl;
    }
}

//-----------------------------------------------------------------------
/// Go through two 1D arrays and compare items item by item
//-----------------------------------------------------------------------

template<class T>
inline void compare_item_by_item(const blitz::Array<T, 1>& A, const blitz::Array<T, 1>& B, double tol, bool show_all=false) {
    compare_shapes(A, B);
    double max = 0;
    for(int i = 0; i < A.extent(blitz::firstDim); i++) {
        double diff = fabs(A(i) - B(i));
        if (diff > max) { max = diff; }
        if( diff >= tol or show_all ) {
            std::cerr << "@[" << i << "] " << A(i) << " != " << B(i) << " -- " << diff << std::endl;
        }
    }
    std::cerr << "Maximum difference = " << max << std::endl;
}

//-----------------------------------------------------------------------
/// Go through two 2D arrays and compare items item by item
//-----------------------------------------------------------------------

template<class T>
inline void compare_item_by_item(const blitz::Array<T, 2>& A, const blitz::Array<T, 2>& B, double tol, bool show_all=false) {
    compare_shapes(A, B);
    double max = 0;
    for(int i = 0; i < A.extent(blitz::firstDim); i++) {
        for(int j = 0; j < A.extent(blitz::secondDim); j++) {
            double diff = fabs(A(i,j) - B(i,j));
            if (diff > max) { max = diff; }
            if( diff >= tol or show_all ) {
                std::cerr << "@[" << i << ", " << j << "] " << A(i,j) << " != " << B(i,j) << " -- " << diff << std::endl;
            }
        }
    }
    std::cerr << "Maximum difference = " << max << std::endl;
}

//-----------------------------------------------------------------------
/// Go through two 3D arrays and compare items item by item
//-----------------------------------------------------------------------

template<class T>
inline void compare_item_by_item(const blitz::Array<T, 3>& A, const blitz::Array<T, 3>& B, double tol, bool show_all=false) {
    compare_shapes(A, B);
    double max = 0;
    for(int i = 0; i < A.extent(blitz::firstDim); i++) {
        for(int j = 0; j < A.extent(blitz::secondDim); j++) {
            for(int k = 0; k < A.extent(blitz::thirdDim); k++) {
                double diff = fabs(A(i,j,k) - B(i,j,k));
                if (diff > max) { max = diff; }
                if( diff >= tol or show_all ) {
                    std::cerr << "@[" << i << ", " << j << ", " << k << "] " << A(i,j,k) << " != " << B(i,j,k) << " -- " << diff << std::endl;
                }
            }
        }
    }
    std::cerr << "Maximum difference = " << max << std::endl;
}

//-----------------------------------------------------------------------
/// Go through two 4D arrays and compare items item by item
//-----------------------------------------------------------------------

template<class T>
inline void compare_item_by_item(const blitz::Array<T, 4>& A, const blitz::Array<T, 4>& B, double tol, bool show_all=false) {
    compare_shapes(A, B);
    double max = 0;
    for(int i = 0; i < A.extent(blitz::firstDim); i++) {
        for(int j = 0; j < A.extent(blitz::secondDim); j++) {
            for(int k = 0; k < A.extent(blitz::thirdDim); k++) {
                for(int m = 0; m < A.extent(blitz::fourthDim); m++) {
                    double diff = fabs(A(i,j,k,m) - B(i,j,k,m));
                    if (diff > max) { max = diff; }
                    if( diff >= tol or show_all ) {
                        std::cerr << "@[" << i << ", " << j << ", " << k << ", " << m << "] " << A(i,j,k,m) << " != " << B(i,j,k,m) << " -- " << diff << std::endl;
                    }
                }
            }
        }
    }
    std::cerr << "Maximum difference = " << max << std::endl;
}

//-----------------------------------------------------------------------
/// Mark as test as "long". We only run these tests when doing
/// "long_check", the normal "check" doesn't run these.  This is
/// handled by checking for the environment variable L2_FP_LONG_CHECK
/// (so when debugging long tests you need to define this environment
/// variable)
//-----------------------------------------------------------------------

#define is_long_test() \
  if(!getenv("L2_FP_LONG_CHECK") && !getenv("L2_FP_REALLY_LONG_CHECK")) { \
    BOOST_WARN_MESSAGE(false, "Skipping long test. To run, make long_check"); \
    return; \
  }

#define is_really_long_test() \
  if(!getenv("L2_FP_REALLY_LONG_CHECK")) { \
    BOOST_WARN_MESSAGE(false, "Skipping really long test. To run, make really_long_check"); \
    return; \
  }

//-----------------------------------------------------------------------
/// Mark as test as "timing". We only run these tests when doing
/// "timing_check", the normal "check" doesn't run these.  This is
/// handled by checking for the environment variable L2_FP_TIMING_CHECK
/// (so when debugging long tests you need to define this environment
/// variable)
//-----------------------------------------------------------------------

#define is_timing_test() \
  if(!getenv("L2_FP_TIMING_CHECK")) { \
    BOOST_WARN_MESSAGE(false, "Skipping timing test. To run, make timing_check"); \
    return; \
  }

#endif
