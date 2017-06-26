// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "aerosol_property_imp_base.h"
#include "sub_state_vector_array.h"
#include "absorber.h"
#include "temperature.h"
#include "altitude.h"
%}

%base_import(aerosol_property)
%base_import(sub_state_vector_array)

%fp_shared_ptr(FullPhysics::AerosolPropertyImpBase);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::AerosolProperty>);
namespace FullPhysics {
%template(SubStateVectorArrayAerosolProperty) 
     FullPhysics::SubStateVectorArray<AerosolProperty>;

// Allow these classes to be derived from in Python.
%feature("director") AerosolPropertyImpBase;

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

class AerosolPropertyImpBase: public SubStateVectorArray<AerosolProperty> {
public:
  // From AerosolPropertyImpBase
  virtual ~AerosolPropertyImpBase();
  virtual boost::shared_ptr<AerosolProperty> clone() const = 0;
  virtual boost::shared_ptr<AerosolProperty> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<RelativeHumidity>& Rh) const = 0;
  virtual ArrayAd<double, 1> extinction_coefficient_each_layer(double wn) 
    const = 0;
  virtual ArrayAd<double, 1> scattering_coefficient_each_layer(double wn) 
    const = 0;
  virtual ArrayAd<double, 3> 
  phase_function_moment_each_layer(double wn, int nmom = -1, 
				   int nscatt = -1) const = 0;
  virtual std::string desc() const;
  %sub_state_virtual_func(AerosolProperty);
  %python_attribute(aerosol_parameter, blitz::Array<double, 1>);
  %python_attribute(aerosol_parameter_uncertainty, blitz::Array<double, 1>);
protected:
  void init(const blitz::Array<double, 1>& Coeff, 
	    const blitz::Array<bool, 1>& Used_flag);
  AerosolPropertyImpBase();
  AerosolPropertyImpBase(const blitz::Array<double, 1>& Coeff, 
			 const blitz::Array<bool, 1>& Used_flag);
};
}

