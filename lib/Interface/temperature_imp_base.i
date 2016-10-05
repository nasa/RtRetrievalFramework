// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"
%{
#include "temperature_imp_base.h"
%}

%base_import(temperature)
%base_import(sub_state_vector_array)

%fp_shared_ptr(FullPhysics::TemperatureImpBase);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::Temperature>);
namespace FullPhysics {
%template(SubStateVectorArrayTemperature) 
     FullPhysics::SubStateVectorArray<Temperature>;

// Allow these classes to be derived from in Python.
%feature("director") TemperatureImpBase;

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

class TemperatureImpBase: public SubStateVectorArray<Temperature> {
public:
  // From PressureImpBase
  virtual ~TemperatureImpBase();
  virtual boost::shared_ptr<Temperature> clone() const = 0;
  virtual boost::shared_ptr<Temperature> 
  clone(const boost::shared_ptr<Pressure>& Press) const = 0;
  virtual AutoDerivativeWithUnit<double> 
  temperature(const AutoDerivativeWithUnit<double>& Press) const;

  %sub_state_virtual_func(Temperature);
protected:
  mutable bool cache_stale;
  mutable boost::function<AutoDerivative<double>(const AutoDerivative<double>&)> tgrid;
  virtual void calc_temperature_grid() const = 0;
  TemperatureImpBase(const blitz::Array<double, 1>& Coeff, 
		     const blitz::Array<bool, 1>& Used_flag,
		     const boost::shared_ptr<Pressure>& Press,
		     bool Mark_according_to_press = true,
	             int Pdep_start = 0);
};
}

