#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H
#include <blitz/array.h>
#include <limits>

namespace FullPhysics {
/** \defgroup LinearAlgebra Linear algebra routines */
/*@{*/

/// \todo Move this to a better location
inline double sqr(double x) {return x * x;}

//-----------------------------------------------------------------------
/// Ensure that a given blitz::Array is contiguous, not reversed, and 
/// in fortran ColumnMajorArray format. If it already is, then we just
/// return the array. Otherwise, we create a copy of this suitable for 
/// passing to Fortran.
//-----------------------------------------------------------------------

template<class T, int D> blitz::Array<T, D> 
to_fortran(const blitz::Array<T, D> &In)
{
  bool is_ok = In.isStorageContiguous();
  for(int i = 0; i < D; ++i)
    if(In.ordering(i) != i ||
       !In.isRankStoredAscending(i))
      is_ok = false;
  if(is_ok)
    return In;
  else {
    blitz::Array<T, D> res(In.shape(), blitz::ColumnMajorArray<D>());
    res = In;
    return res;
  }
}

//-----------------------------------------------------------------------
/// Ensure that a given blitz::Array is contiguous, not reversed, and 
/// in C RowMajorArray format. If it already is, then we just
/// return the array. Otherwise, we create a copy of this in C order.
//-----------------------------------------------------------------------

template<class T, int D> blitz::Array<T, D> 
to_c_order(const blitz::Array<T, D> &In)
{
  bool is_ok = In.isStorageContiguous();
  for(int i = 0; i < D; ++i)
    if(In.ordering(i) != D - i - 1||
       !In.isRankStoredAscending(i))
      is_ok = false;
  if(is_ok)
    return In;
  else {
    blitz::Array<T, D> res(In.shape());
    res = In;
    return res;
  }
}

//-----------------------------------------------------------------------
/// Ensure that a given blitz::Array is contiguous, not reversed, and 
/// in C RowMajorArray format. If it already is, then we just
/// return the array. Otherwise, we create a copy of this in C order.
///
/// Unfortunately, this removes the "const" from the calling argument,
/// so we have a version of this where you explicitly say this ok.
//-----------------------------------------------------------------------

template<class T, int D> blitz::Array<T, D> 
to_c_order_const(const blitz::Array<T, D> &In)
{
  bool is_ok = In.isStorageContiguous();
  for(int i = 0; i < D; ++i)
    if(In.ordering(i) != D - i - 1||
       !In.isRankStoredAscending(i))
      is_ok = false;
  if(is_ok)
    return In;
  else {
    blitz::Array<T, D> res(In.shape());
    res = In;
    return res;
  }
}

//-----------------------------------------------------------------------
/// If matrix is in Fortran order, resize it. Otherwise, set it to a
/// Fortran matrix of the given size.
//-----------------------------------------------------------------------

template<class T> void fortran_resize(blitz::Array<T, 1>& M, int sz1)
{
  const int D = 1;
  bool is_ok = M.isStorageContiguous();
  for(int i = 0; i < D; ++i)
    if(M.ordering(i) != i ||
       !M.isRankStoredAscending(i))
      is_ok = false;
  if(is_ok)
    M.resize(sz1);
  else
    M.reference(blitz::Array<T, D>(sz1, blitz::ColumnMajorArray<D>()));
}

//-----------------------------------------------------------------------
/// If matrix is in Fortran order, resize it. Otherwise, set it to a
/// Fortran matrix of the given size.
//-----------------------------------------------------------------------

template<class T> void fortran_resize(blitz::Array<T, 2>& M, int sz1, int sz2)
{
  const int D = 2;
  bool is_ok = M.isStorageContiguous();
  for(int i = 0; i < D; ++i)
    if(M.ordering(i) != i ||
       !M.isRankStoredAscending(i))
      is_ok = false;
  if(is_ok)
    M.resize(sz1, sz2);
  else
    M.reference(blitz::Array<T, D>(sz1, sz2, blitz::ColumnMajorArray<D>()));
}

//-----------------------------------------------------------------------
/// If matrix is in Fortran order, resize it. Otherwise, set it to a
/// Fortran matrix of the given size.
//-----------------------------------------------------------------------

template<class T> void fortran_resize(blitz::Array<T, 3>& M, int sz1, int sz2, 
                                      int sz3)
{
  const int D = 3;
  bool is_ok = M.isStorageContiguous();
  for(int i = 0; i < D; ++i)
    if(M.ordering(i) != i ||
       !M.isRankStoredAscending(i))
      is_ok = false;
  if(is_ok)
    M.resize(sz1, sz2, sz3);
  else
    M.reference(blitz::Array<T, D>(sz1, sz2, sz3, 
                                   blitz::ColumnMajorArray<D>()));
}

//-----------------------------------------------------------------------
/// Ensure that a given blitz::Array is contiguous, not reversed, and 
/// in fortran ColumnMajorArray format. If it already is, then we just
/// return the array. Otherwise, we create a copy of this suitable for 
/// passing to Fortran.
///
/// Unfortunately, this removes the "const" from the calling argument,
/// so we have a version of this where you explicitly say this ok.
//-----------------------------------------------------------------------

template<class T, int D> blitz::Array<T, D> 
to_fortran_const(const blitz::Array<T, D> &In)
{
  bool is_ok = In.isStorageContiguous();
  for(int i = 0; i < D; ++i)
    if(In.ordering(i) != i ||
       !In.isRankStoredAscending(i))
      is_ok = false;
  if(is_ok)
    return In;
  else {
    blitz::Array<T, D> res(In.shape(), blitz::ColumnMajorArray<D>());
    res = In;
    return res;
  }
}

blitz::Array<double, 1> solve_least_squares(const blitz::Array<double, 2>& A, 
                                            const blitz::Array<double, 1>& B, 
                                            double Rcond = 1e-12);

blitz::Array<double, 1> solve_least_squares_qr(const blitz::Array<double, 2>& A, 
                                               const blitz::Array<double, 1>& B);

void svd(const blitz::Array<double, 2>& A, blitz::Array<double, 1>& S,
         blitz::Array<double, 2>& U, blitz::Array<double, 2>& VT);

blitz::Array<double, 2> cholesky_decomposition(const blitz::Array<double, 2>& A);

blitz::Array<double, 2> generalized_inverse(const blitz::Array<double, 2>& A, 
                    double Rcond = std::numeric_limits<double>::epsilon());
blitz::Array<double, 2> generalized_inverse(const blitz::Array<double, 2>& A, 
                    const blitz::Array<bool, 1>& Zero_unused,
                    double Rcond = std::numeric_limits<double>::epsilon());


blitz::Array<double, 2> multifit_covar(const blitz::Array<double, 2>& J, 
                                       double Eps_rel = std::numeric_limits<double>::epsilon());
bool multifit_test_gradient(const blitz::Array<double, 1>& g, 
                            double Eps_abs = std::numeric_limits<double>::epsilon());
bool multifit_test_delta(const blitz::Array<double, 1>& dx, const blitz::Array<double, 1>& x,
                         double Eps_abs = std::numeric_limits<double>::epsilon(),
                         double Eps_rel = std::numeric_limits<double>::epsilon());

/*@}*/
}

#endif
