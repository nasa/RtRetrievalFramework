// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"
%{
#include "sub_state_vector_array.h"
#include "absorber.h"
#include "temperature.h"
#include "altitude.h"
%}
%base_import(state_vector)
%base_import(observer)
%import "sub_state_vector_array.i"
%import "pressure.i"
%import "temperature.i"
%import "altitude.i"
%import "absorber_vmr.i"

%fp_shared_ptr(FullPhysics::Absorber)

namespace FullPhysics {
  class Absorber;
}

%nodefaultctor FullPhysics::SubStateVectorArray<FullPhysics::Absorber>;
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Absorber>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Absorber>)
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::Absorber>);

namespace FullPhysics {
%template(ObservableAbsorber) FullPhysics::Observable<Absorber>;
%template(ObserverAbsorber) FullPhysics::Observer<Absorber>;

// Allow these classes to be derived from in Python.
%feature("director") Absorber;

// Note, a class that is derived from in python needs to declare every virtual function that
// can be called on it, even if all that happens is the base class
// to a director is called. This is because this class is used to
// create the SwigDirector class, and this class needs each of the member functions to
// direct things properly. It is *not* necessary to add these function to the underlying
// C++, only that you declare them here.

class Absorber : virtual public StateVectorObserver, 
		 public Observable<Absorber> {
public:
  virtual ~Absorber();
  virtual void add_observer(Observer<Absorber>& Obs);
  virtual void remove_observer(Observer<Absorber>& Obs);
  std::string print_to_string() const;
  %python_attribute(number_species, virtual int);
  virtual std::string gas_name(int Species_index) const = 0;
  virtual int gas_index(const std::string& Name) const;
  
  virtual ArrayAd<double, 2> 
  optical_depth_each_layer(double wn, int spec_index) const = 0;
  virtual AutoDerivative<double> xgas(const std::string& Gas_name) const = 0;
  virtual boost::shared_ptr<AbsorberVmr> absorber_vmr(const std::string& gas_name) const = 0;
  virtual boost::shared_ptr<Absorber> clone() const = 0;
  virtual boost::shared_ptr<Absorber> clone
    (const boost::shared_ptr<Pressure>& Press,
     const boost::shared_ptr<Temperature>& Temp,
     const std::vector<boost::shared_ptr<Altitude> >& Alt) const = 0;
  virtual void print(std::ostream& Os) const;

  // Functions so StateVectorObserver base class, see note above
  virtual void notify_update(const StateVector& Sv);
  virtual void mark_used(const StateVector& Sv, 
			 blitz::Array<bool, 1>& Used) const;
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;
  virtual void notify_add(StateVector& Observed_object);
  virtual void notify_remove(StateVector& Observed_object);
};
} // Rename below doesn't like being in a namespace
%feature("director") FullPhysics::SubStateVectorArray<FullPhysics::Absorber>;

// Note, a class that is derived from in python needs to declare every
// virtual function that can be called on it, even if all that happens
// is the base class. This is because this class is used to create the
// SwigDirector class, and this class needs each of the member
// functions to direct things properly. It is *not* necessary to add
// these function to the underlying C++, only that you declare them
// here.
//
// Since We need to add some member functions to this template, so we
// need to explictly expand it in SWIG. This takes two parts, the
// rename here and the declaration below.
%rename(SubStateVectorAbsorber) FullPhysics::SubStateVectorArray<FullPhysics::Absorber>;

namespace FullPhysics {

class FullPhysics::SubStateVectorArray<FullPhysics::Absorber>: 
    public Absorber,
    public SubStateVectorObserver {
public:
  SubStateVectorArray<FullPhysics::Absorber>(const blitz::Array<double, 1>& Coeff, 
		      const blitz::Array<bool, 1>& Used_flag);
  SubStateVectorArray<FullPhysics::Absorber>(const blitz::Array<double, 1>& Coeff, 
		      const blitz::Array<bool, 1>& Used_flag,
		      const boost::shared_ptr<Pressure>& Press);
  SubStateVectorArray<FullPhysics::Absorber>();

  // Functions from Absorber base class, see note above.
  virtual int number_species() const;
  virtual std::string gas_name(int Species_index) const = 0;
  virtual int gas_index(const std::string& Name) const;
  
  virtual ArrayAd<double, 2> 
  optical_depth_each_layer(double wn, int spec_index) const = 0;
  virtual void print(std::ostream& Os) const;
  virtual void add_observer(Observer<Absorber>& Obs); 
  virtual void remove_observer(Observer<Absorber>& Obs);
  virtual boost::shared_ptr<Absorber> clone() const = 0;
  virtual boost::shared_ptr<Absorber> clone
    (const boost::shared_ptr<Pressure>& Press,
     const boost::shared_ptr<Temperature>& Temp,
     const std::vector<boost::shared_ptr<Altitude> >& Alt) const = 0;
  virtual void notify_update(const StateVector& Sv);
  virtual void mark_used(const StateVector& Sv, 
			 blitz::Array<bool, 1>& Used) const;
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;
  virtual void notify_add(StateVector& Observed_object);
  virtual void notify_remove(StateVector& Observed_object);

  void init(const blitz::Array<double, 1>& Coeff, 
	    const blitz::Array<bool, 1>& Used_flag,
	    const boost::shared_ptr<Pressure>& Press =
	    boost::shared_ptr<Pressure>());
  void mark_used_sub(blitz::Array<bool, 1>& Used) const;
  virtual std::string state_vector_name_i(int i) const;
  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) 
    const;
  
  virtual void update_sub_state(const ArrayAd<double, 1>& Sv_sub,
				const blitz::Array<double, 2>& Cov);
  virtual void update_sub_state_hook();
  ArrayAd<double, 1> coefficient() const;
  blitz::Array<bool, 1> used_flag_value() const;
  boost::shared_ptr<Pressure> pressure() const;
};

}


