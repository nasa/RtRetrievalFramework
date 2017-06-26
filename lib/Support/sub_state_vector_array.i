// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "sub_state_vector_array.h"
#include "pressure.h"
%}

%import "state_vector.i"
%import "pressure.i"

namespace FullPhysics {
class Pressure;
template<class Base> class SubStateVectorArray: 
    public Base,
    public SubStateVectorObserver {
public:
  SubStateVectorArray(const blitz::Array<double, 1>& Coeff, 
		      const blitz::Array<bool, 1>& Used_flag);
  SubStateVectorArray(const blitz::Array<double, 1>& Coeff, 
		      const blitz::Array<bool, 1>& Used_flag,
		      const boost::shared_ptr<Pressure>& Press,
		      bool Mark_according_to_press = true,
		      int Pdep_start = 0);
  SubStateVectorArray();
  void init(const blitz::Array<double, 1>& Coeff, 
	    const blitz::Array<bool, 1>& Used_flag);
  void init(const blitz::Array<double, 1>& Coeff, 
	    const blitz::Array<bool, 1>& Used_flag,
	    const boost::shared_ptr<Pressure>& Press,
	    bool Mark_according_to_press = true,
	    int Pdep_start = 0);
  virtual ~SubStateVectorArray();
  void mark_used_sub(blitz::Array<bool, 1>& Used) const;
  virtual std::string state_vector_name_i(int i) const;
  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) 
    const;
  
  virtual void update_sub_state(const ArrayAd<double, 1>& Sv_sub,
				const blitz::Array<double, 2>& Cov);
  virtual void update_sub_state_hook();
  %python_attribute(coefficient, ArrayAd<double, 1>);
  %python_attribute(used_flag_value, blitz::Array<bool, 1>);
  %python_attribute(statevector_covariance, blitz::Array<double, 2>);
  %python_attribute(pressure, boost::shared_ptr<Pressure>);
protected:
  ArrayAd<double, 1> coeff;
  boost::shared_ptr<Pressure> press;
  blitz::Array<bool, 1> used_flag;
  blitz::Array<double, 2> cov;
  bool mark_according_to_press;
  int pdep_start;
};
}

// When we create directors, at least for SWIG 2.0.4 we need to explicitly
// list every virtual function. We have lots of classes that derive from
// SubStateVectorArray, so we have this utility macro to define all those
// virtual functions.
%define %sub_state_virtual_func(TYPE)
virtual void add_observer(Observer<TYPE>& Obs);
virtual void remove_observer(Observer<TYPE>& Obs);
virtual void update_sub_state_hook();
virtual void print(std::ostream& Os) const;
%python_attribute(desc, std::string)
virtual void mark_used(const StateVector& Sv, 
		 blitz::Array<bool, 1>& Used) const;
virtual void state_vector_name(const StateVector& Sv, 
				 blitz::Array<std::string, 1>& Sv_name) const;
virtual void notify_update(const StateVector& Observed_object);
virtual void notify_add(StateVector& Observed_object);
virtual void notify_remove(StateVector& Observed_object);
virtual void update_sub_state(const ArrayAd<double, 1>& Sv_sub,
				const blitz::Array<double, 2>& Cov_sub);
virtual std::string state_vector_name_i(int i) const;
virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) 
    const;
%enddef
