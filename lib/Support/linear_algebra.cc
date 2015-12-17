#include "linear_algebra.h"
#include "fp_gsl_matrix.h"
#include "fp_exception.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_multifit_nlin.h"

using namespace FullPhysics;
using namespace blitz;

/** \addtogroup LinearAlgebra */
/*@{*/

//-----------------------------------------------------------------------
/// This solves the least squares system A*x = b, returning the
/// minimum norm solution. It is perfectly ok for A to be rank
/// deficient.
///
/// \param A a M x N matrix, which may be rank deficient.
/// \param B the right hand side, a M vector
/// \param Rcond a singular value <= Rcond * max singular value is
///     treated as 0. You can pass Rcond < 0 to use machine precision
///     if desired.
/// \return The vector X, which is a N vector
//-----------------------------------------------------------------------

Array<double, 1> 
FullPhysics::solve_least_squares(const blitz::Array<double, 2>& A, 
                                 const blitz::Array<double, 1>& B, 
                                 double Rcond)
{
  int n = A.cols();
  Array<double, 1> s;
  Array<double, 2> u;
  Array<double, 2> v;
  svd(A, s, u, v);
  v.transposeSelf(secondDim, firstDim);
  // Set very small values to zero.
  for(int i = 0; i < s.rows(); ++i)
    if(s(i) <= s(0) * Rcond)
      s(i) = 0.0;
  GslVector sgsl(s);
  GslMatrix ugsl(u);
  GslMatrix vgsl(v);
  GslVector bgsl(const_cast<Array<double, 1>&>(B));
  Array<double, 1> x(n);
  GslVector xgsl(x);
  int status = gsl_linalg_SV_solve(ugsl.gsl(), vgsl.gsl(), sgsl.gsl(), 
                                   bgsl.gsl(), xgsl.gsl());
  gsl_check(status);
  return xgsl.blitz_array();
}

//-----------------------------------------------------------------------
/// This finds the least squares solution to the overdetermined 
/// system A x = b where the matrix A has more rows than columns. 
/// The least squares solution minimizes the Euclidean norm of the 
/// residual, ||Ax - b||.
/// The solution is determined using a QR decomposition of A.
///
/// \param A a M x N matrix, which may be rank deficient.
/// \param B the right hand side, a M vector
/// \return The vector X, which is a N vector
//-----------------------------------------------------------------------

Array<double, 1> 
FullPhysics::solve_least_squares_qr(const blitz::Array<double, 2>& A, 
                                    const blitz::Array<double, 1>& B)
{
  GslMatrix Agsl(const_cast<Array<double, 2>&>(A));
  GslVector Bgsl(const_cast<Array<double, 1>&>(B));

  // tau must be the min(M, N)
  Array<double, 1> tau;
  tau.resize(min(A.rows(), A.cols()));
  GslVector taugsl(tau);

  int status = gsl_linalg_QR_decomp (Agsl.gsl(), taugsl.gsl());
  gsl_check(status);
 
  Array<double, 1> X(A.cols());
  GslVector Xgsl(X);

  Array<double, 1> residual(A.rows());
  GslVector residualgsl(residual);

  status = gsl_linalg_QR_lssolve (Agsl.gsl(), taugsl.gsl(), Bgsl.gsl(), Xgsl.gsl(), residualgsl.gsl());
  gsl_check(status);

  return Xgsl.blitz_array();
}

//-----------------------------------------------------------------------
/// Compute the SVD decomposition of a matrix A = U * SIGMA * V^T
///
/// \param A Matrix to get decomposition for.
/// \param S On exit, Singular values
/// \param U On exit, U matrix. Note that this is the "thin" version
///        of U, this is M x N rather than M x M. The "full" version
///        would just have extra columns of zeros.
/// \param VT On exit, V^T matrix
//-----------------------------------------------------------------------

void FullPhysics::svd(const blitz::Array<double, 2>& A, 
                      blitz::Array<double, 1>& S,
                      blitz::Array<double, 2>& U, 
                      blitz::Array<double, 2>& VT)
{
  int m = A.rows();
  int n = A.cols();
  range_min_check(m, n);
  S.resize(n);
  U.reference(Array<double,2>(A.shape()));
  U = A;                        // This is calculated in place.
  VT.reference(Array<double,2>(shape(n,n)));
  Array<double, 1> work(n);
  GslMatrix ugsl(U);
  GslVector sgsl(S);
  GslMatrix vtgsl(VT);
  GslVector workgsl(work);
  int status = gsl_linalg_SV_decomp(ugsl.gsl(), vtgsl.gsl(), sgsl.gsl(), 
                                    workgsl.gsl());
  gsl_check(status);
  VT.transposeSelf(secondDim, firstDim);
}

//-----------------------------------------------------------------------
/// This returns the generalized inverse of the given matrix. This is
/// the same as the inverse, except that if A is rank deficient we set
/// the SVD value to 0.
///
/// \param A Matrix to get inverse of 
/// \param Rcond We set singular values < Rcond * max singular value
/// to 0.
//-----------------------------------------------------------------------

Array<double, 2> FullPhysics::generalized_inverse(const 
          blitz::Array<double, 2>& A, double Rcond)
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  blitz::Array<double, 1> s;
  blitz::Array<double, 2> u, vt;
  svd(A, s, u, vt);
  Array<double, 2> res(vt.cols(), u.rows());
  res = sum(where(s(i3) > s(0) * Rcond, vt(i3, i1) * u(i2, i3) / s(i3), 0), 
            i3);
  return res;
}

//-----------------------------------------------------------------------
/// This returns the generalized inverse of the given matrix. This is
/// the same as the inverse, except that if A is rank deficient we set
/// the SVD value to 0.
///
/// This variation takes a boolean array indicating which parameters
/// are to be held at zero.
///
/// \param A Matrix to get inverse of 
/// \param Zero_unused Boolean array indicating parameters that are
///   unused and should be held to zero
/// \param Rcond We set singular values < Rcond * max singular value
///     to 0.
//-----------------------------------------------------------------------

Array<double, 2> FullPhysics::generalized_inverse(const 
                  blitz::Array<double, 2>& A, const blitz::Array<bool, 1>& Zero_unused, 
                  double Rcond)
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  if(A.rows() != A.cols())
    throw Exception("This version of generalized_inverse can only be used with square matrix");
  blitz::Array<double, 1> s;
  blitz::Array<double, 2> u, vt;
  int num_zero = count(Zero_unused);
  blitz::Array<double, 2> asub(A.rows() - num_zero, A.cols() - num_zero);
  int isub = 0;
  for(int i = 0; i < A.rows(); ++i) 
    if(!Zero_unused(i)) {
      int jsub = 0;
      for(int j = 0; j < A.cols(); ++j)
        if(!Zero_unused(j)) {
          asub(isub, jsub) = A(i, j);
          ++jsub;
        }
      ++isub;
    }
  svd(asub, s, u, vt);
  Array<double, 2> ressub(vt.cols(), u.rows());
  ressub = sum(where(s(i3) > s(0) * Rcond, vt(i3, i1) * u(i2, i3) / s(i3), 0), 
               i3);
  Array<double, 2> res(A.shape());
  res = 0;
  isub = 0;
  for(int i = 0; i < A.rows(); ++i) 
    if(!Zero_unused(i)) {
      int jsub = 0;
      for(int j = 0; j < A.cols(); ++j)
        if(!Zero_unused(j)) {
          res(i,j) = ressub(isub, jsub);
          ++jsub;
        }
      ++isub;
    }
  return res;
}


//-----------------------------------------------------------------------
/// This calculates the Cholesky Decompostion of A so that 
/// A = L L^T. This returns L. A must be symmetric positive definite.
//-----------------------------------------------------------------------

blitz::Array<double, 2> FullPhysics::cholesky_decomposition(const blitz::Array<double, 2>& A)
{
  Array<double, 2> t(A.shape());
  GslMatrix rgsl(t);
  Array<double, 2> res(rgsl.blitz_array());
  res = A;
  int status = gsl_linalg_cholesky_decomp(rgsl.gsl());
  gsl_check(status);
  for(int i = 0; i < res.rows(); ++i)
    for(int j = i + 1; j < res.cols(); ++j)
      res(i, j) = 0;
  return res;
}

blitz::Array<double, 2> FullPhysics::multifit_covar(const blitz::Array<double, 2>& J, double Eps_rel)
{
  Array<double, 2> covar(J.cols(), J.cols());
  int status = gsl_multifit_covar(GslMatrix(const_cast< Array<double, 2>& >(J)).gsl(), Eps_rel, GslMatrix(covar).gsl());
  gsl_check(status);
  return covar;
}

bool FullPhysics::multifit_test_gradient(const blitz::Array<double, 1>& g, double Eps_abs)
{
  return gsl_multifit_test_gradient(GslVector(const_cast< Array<double, 1>& >(g)).gsl(), Eps_abs) == GSL_SUCCESS;
}

bool FullPhysics::multifit_test_delta(
    const blitz::Array<double, 1>& dx, const blitz::Array<double, 1>& x, double Eps_abs, double Eps_rel)
{
  return gsl_multifit_test_delta( GslVector(const_cast< Array<double, 1>& >(dx)).gsl(), 
                                  GslVector(const_cast< Array<double, 1>& >(x)).gsl(), 
                                  Eps_abs, Eps_rel) == GSL_SUCCESS;
}

/*@}*/
