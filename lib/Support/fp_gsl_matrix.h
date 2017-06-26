#ifndef FP_GSL_MATRIX_H
#define FP_GSL_MATRIX_H
#include <gsl/gsl_matrix.h>
#include <blitz/array.h>
#include <boost/utility.hpp>

namespace FullPhysics {

/****************************************************************//**
  This provides thin wrapper around the GNU Scientific Library
  gsl_matrix.

  The GSL is a pretty complete scientific library. However, it is C based
  and out of the box doesn't place nicely with other C++ classes (such
  as blitz::Array). We provide some thin wrapper code around the GSL
  code to interface these two classes. This provides both a
  blitz::Array and a gsl_matrix view of the same underlying data.
*******************************************************************/

class GslMatrix : public boost::noncopyable {
public:
//-----------------------------------------------------------------------
/// Default constructor.
//-----------------------------------------------------------------------

  GslMatrix() : gsl_matrix_(0) {}

//-----------------------------------------------------------------------
/// Use data owned by gsl_matrix. This data can either have ownership
/// passed to this class (in which case we delete it when done with
/// it), or just a reference (in which case the lifetime is handled
/// outside of this class).
//-----------------------------------------------------------------------

  GslMatrix(gsl_matrix* M, bool Owned = true) : gsl_matrix_(0) 
  {reset(M, Owned);}

//-----------------------------------------------------------------------
/// Use data owned by blitz::Array. Note that if this data
/// isStorageContiguous(), and in C order then we use the data in
/// place. If it isn't 
/// be make a contiguous copy of it. This means that in one case the
/// original matrix will be modified if the gsl_matrix view of the
/// data is modified, and the second case it won't. If you don't want
/// this behaviour, you can form a copy before you pass it to this
/// constructor. 
//-----------------------------------------------------------------------

  GslMatrix(blitz::Array<double, 2>& M) : gsl_matrix_(0) {reset(M);}
  ~GslMatrix();
  void reset(gsl_matrix* M, bool Owned = true);
  void reset(blitz::Array<double, 2>& M);

//-----------------------------------------------------------------------
/// Return gsl_matrix look at data.
//-----------------------------------------------------------------------

  const gsl_matrix* gsl() const { return gsl_matrix_; }

//-----------------------------------------------------------------------
/// Return gsl_matrix look at data.
//-----------------------------------------------------------------------

  gsl_matrix* gsl() { return gsl_matrix_; }

//-----------------------------------------------------------------------
/// Return blitz::Array look at data.
//-----------------------------------------------------------------------

  const blitz::Array<double, 2>& blitz_array() const {return blitz_array_;}

//-----------------------------------------------------------------------
/// Return blitz::Array look at data.
//-----------------------------------------------------------------------

  blitz::Array<double, 2>& blitz_array() {return blitz_array_;}
public:
  blitz::Array<double, 2> blitz_array_;
  gsl_matrix* gsl_matrix_;
  bool owned_;			// This only applies to gsl_matrix_,
				// blitz_array manages itself.
  bool block_owned_;
};

/****************************************************************//**
  This provides thin wrapper around the GNU Scientific Library
  gsl_vector.

  The GSL is a pretty complete scientific library. However, is C based
  and out of the box doesn't place nicely with other C++ classes (such
  as blitz::Array). We provide some thin wrapper code around the GSL
  code to interface these two classes. This provides both a
  blitz::Array and a gsl_vector view of the same underlying data.
*******************************************************************/

class GslVector : public boost::noncopyable {
public:
//-----------------------------------------------------------------------
/// Default constructor.
//-----------------------------------------------------------------------

  GslVector() : gsl_vector_(0) {}

//-----------------------------------------------------------------------
/// Use data owned by gsl_vector. This data can either have ownership
/// passed to this class (in which case we delete it when done with
/// it), or just a reference (in which case the lifetime is handled
/// outside of this class).
//-----------------------------------------------------------------------

  GslVector(gsl_vector* M, bool Owned = true) :gsl_vector_(0) 
  { reset(M, Owned);}

//-----------------------------------------------------------------------
/// Use data owned by blitz::Array. Note that if this data
/// isStorageContiguous(), and in C order then we use the data in
/// place. If it isn't 
/// be make a contiguous copy of it. This means that in one case the
/// original vector will be modified if the gsl_vector view of the
/// data is modified, and the second case it won't. If you don't want
/// this behaviour, you can form a copy before you pass it to this
/// constructor. 
//-----------------------------------------------------------------------

  GslVector(blitz::Array<double, 1>& M) 
    :gsl_vector_(0) { reset(M);}
  ~GslVector();
  void reset(gsl_vector* M, bool Owned = true);
  void reset(blitz::Array<double, 1>& M);

//-----------------------------------------------------------------------
/// Return gsl_vector look at data.
//-----------------------------------------------------------------------

  const gsl_vector* gsl() const { return gsl_vector_; }

//-----------------------------------------------------------------------
/// Return gsl_vector look at data.
//-----------------------------------------------------------------------

  gsl_vector* gsl() { return gsl_vector_; }

//-----------------------------------------------------------------------
/// Return blitz::Array look at data.
//-----------------------------------------------------------------------

  const blitz::Array<double, 1>& blitz_array() const {return blitz_array_;}

//-----------------------------------------------------------------------
/// Return blitz::Array look at data.
//-----------------------------------------------------------------------

  blitz::Array<double, 1>& blitz_array() {return blitz_array_;}
public:
  blitz::Array<double, 1> blitz_array_;
  gsl_vector* gsl_vector_;
  bool owned_;			// This only applies to gsl_vector_,
				// blitz_array manages itself.
  bool block_owned_;
};
}
#endif
