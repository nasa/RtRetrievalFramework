// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "stokes_coefficient_imp_base.h"
%}

%base_import(stokes_coefficient)
%import "sub_state_vector_array.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::StokesCoefficientImpBase);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::StokesCoefficient>)

namespace FullPhysics {
%template(SubStateVectorArrayStokesCoefficient) 
     FullPhysics::SubStateVectorArray<StokesCoefficient>;

// Note, at least for SWIG 2.0.4 a class that is derived from in python 
// needs to declare every virtual function that can be called on it, even 
// if all that happens is the base class to a director is called. This is 
// because this class is used to create the SwigDirector class, and this 
// class needs each of the member functions to direct things properly. It 
// is *not* necessary to add these function to the underlying
// C++, only that you declare them here.
//
// If you miss something, then you will get something like a recursion
// error in python when a virtual function is used that isn't explicitly
// listed here.
//
// This seems like a bug in 2.0.4, if SWIG needs all the member functions
// it should know to make them itself. So perhaps a future version of SWIG
// won't have this same constraint. But for now, this is required.

class StokesCoefficientImpBase: public SubStateVectorArray<StokesCoefficient> {
public:
  // From StokesCoefficientImpBase
  virtual ~StokesCoefficientImpBase();
  virtual boost::shared_ptr<StokesCoefficient> clone() const = 0;
  %python_attribute_derived(stokes_coefficient, ArrayAd<double, 2>);
  %sub_state_virtual_func(StokesCoefficient);
protected:
  mutable bool cache_stale;
  mutable ArrayAd<double, 2> stokes_coeff;
  virtual void calc_stokes_coeff() const = 0;
  StokesCoefficientImpBase();
  StokesCoefficientImpBase(const blitz::Array<double, 1>& Coeff, 
			   const blitz::Array<bool, 1>& Used_flag);
};
}

