#include "fp_gsl_matrix.h"
#include "linear_algebra.h"
#include <cstdlib>

using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Reset class to point to new matrix. This data can either have ownership
/// passed to this class (in which case we delete it when done with
/// it), or just a reference (in which case the lifetime is handled
/// outside of this class).
//-----------------------------------------------------------------------

void GslMatrix::reset(gsl_matrix* M, bool Owned)
{
  if(gsl_matrix_ && block_owned_)
    free(gsl_matrix_->block);
  if(gsl_matrix_ && owned_)
    gsl_matrix_free(gsl_matrix_);
  gsl_matrix_ = M;
  blitz::Array<double, 2> a(M->data, blitz::shape(M->size1, M->size2), 
			    blitz::shape(M->tda, 1), blitz::neverDeleteData);
  blitz_array_.reference(a);
  owned_ = Owned;
  block_owned_ = false;
}

//-----------------------------------------------------------------------
/// Reset to Use data owned by blitz::Array. Note that if this data
/// isStorageContiguous(), and in C order then we use the data in
/// place. If it isn't 
/// be make a contiguous copy of it. This means that in one case the
/// original matrix will be modified if the gsl_matrix view of the
/// data is modified, and the second case it won't. If you don't want
/// this behaviour, you can form a copy before you pass it to this
/// constructor. 
//-----------------------------------------------------------------------

void GslMatrix::reset(blitz::Array<double, 2>& M)
{
  if(gsl_matrix_ && block_owned_)
    free(gsl_matrix_->block);
  if(gsl_matrix_ && owned_)
    gsl_matrix_free(gsl_matrix_);
  blitz_array_.reference(to_c_order(M));
  gsl_matrix_ = (gsl_matrix*) malloc(sizeof(gsl_matrix));
  gsl_matrix_->size1 = blitz_array_.rows();
  gsl_matrix_->size2 = blitz_array_.cols();
  gsl_matrix_->tda = blitz_array_.cols();
  gsl_matrix_->data = blitz_array_.data();
  gsl_matrix_->owner = 0;
  gsl_matrix_->block = (gsl_block*) malloc(sizeof(gsl_block));
  gsl_matrix_->block->size = blitz_array_.size();
  gsl_matrix_->block->data = blitz_array_.data();
  owned_ = true;
  block_owned_ = true;
}

//-----------------------------------------------------------------------
/// Destructor.
//-----------------------------------------------------------------------

GslMatrix::~GslMatrix()
{
  if(gsl_matrix_ && block_owned_)
    free(gsl_matrix_->block);
  if(gsl_matrix_ && owned_)
    gsl_matrix_free(gsl_matrix_);
}

//-----------------------------------------------------------------------
/// Reset class to point to new vector. This data can either have ownership
/// passed to this class (in which case we delete it when done with
/// it), or just a reference (in which case the lifetime is handled
/// outside of this class).
//-----------------------------------------------------------------------

void GslVector::reset(gsl_vector* M, bool Owned)
{
  if(gsl_vector_ && block_owned_)
    free(gsl_vector_->block);
  if(gsl_vector_ && owned_)
    gsl_vector_free(gsl_vector_);
  gsl_vector_ = M;
  blitz::Array<double, 1> a(M->data, blitz::shape(M->size), 
			    blitz::shape(M->stride), blitz::neverDeleteData);
  blitz_array_.reference(a);
  owned_ = Owned;
  block_owned_ = false;
}

//-----------------------------------------------------------------------
/// Reset to use data owned by blitz::Array. Note that if this data
/// isStorageContiguous(), then we use the data in place. If it isn't
/// be make a contiguous copy of it. This means that in one case the
/// original vector will be modified if the gsl_vector view of the
/// data is modified, and the second case it won't. If you don't want
/// this behaviour, you can form a copy before you pass it to this
/// constructor. 
//-----------------------------------------------------------------------

void GslVector::reset(blitz::Array<double, 1>& M)
{
  if(gsl_vector_ && block_owned_)
    free(gsl_vector_->block);
  if(gsl_vector_ && owned_)
    gsl_vector_free(gsl_vector_);
  blitz_array_.reference(to_c_order(M));
  gsl_vector_ = (gsl_vector*) malloc(sizeof(gsl_vector));
  gsl_vector_->size = blitz_array_.rows();
  gsl_vector_->stride = 1;
  gsl_vector_->data = blitz_array_.data();
  gsl_vector_->owner = 0;
  gsl_vector_->block = (gsl_block*) malloc(sizeof(gsl_block));
  gsl_vector_->block->size = blitz_array_.size();
  gsl_vector_->block->data = blitz_array_.data();
  owned_ = true;
  block_owned_ = true;
}

//-----------------------------------------------------------------------
/// Destructor.
//-----------------------------------------------------------------------

GslVector::~GslVector()
{
  if(gsl_vector_ && block_owned_)
    free(gsl_vector_->block);
  if(gsl_vector_ && owned_) 
    gsl_vector_free(gsl_vector_);
}

