// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "radiative_transfer_imp_base.h"
%}
%base_import(radiative_transfer_retrievable)
%base_import(sub_state_vector_array)
%import "spectral_domain.i"
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::RadiativeTransferImpBase);
// Hide the fact that for Python implementations
// we are returning a shared pointer, not the actual object itself
%rename(reflectance) reflectance_ptr;

namespace FullPhysics {

// Allow these classes to be derived from in Python.
%feature("director") RadiativeTransferImpBase;
%feature("nodirector") RadiativeTransferImpBase::reflectance;

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

class RadiativeTransferImpBase: public SubStateVectorArray<RadiativeTransferRetrievable> {
public:
  // From RadiativeTransferImpBase
  virtual ~RadiativeTransferImpBase();
  virtual boost::shared_ptr<RadiativeTransferRetrievable> clone() const = 0;

  // From RadiativeTransfer
  %python_attribute_abstract(number_stokes, int);
  %python_attribute_abstract(number_spectrometer, int);

  // Return a shared_ptr for Spectrum class because for some ungodly reason
  // SWIG doesn't like returning naked Spectrum classes
  virtual boost::shared_ptr<FullPhysics::Spectrum> reflectance_ptr(const FullPhysics::SpectralDomain& Spec_domain, int Spec_index, 
						   bool Skip_jacobian = false) const = 0;

  virtual blitz::Array<double, 2> stokes(const FullPhysics::SpectralDomain& Spec_domain,
					 int Spec_index) const = 0;
  virtual ArrayAd<double, 2> stokes_and_jacobian(const FullPhysics::SpectralDomain& Spec_domain, 
						 int Spec_index) const = 0;

  %sub_state_virtual_func(RadiativeTransferRetrievable);
protected:
  RadiativeTransferImpBase();
  RadiativeTransferImpBase(const blitz::Array<double, 1>& Coeff, 
			   const blitz::Array<bool, 1>& Used_flag);
};
}

