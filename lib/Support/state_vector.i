// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "state_vector.h"
%}
%base_import(generic_object)
%import "observer.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::StateVector);
%fp_shared_ptr(FullPhysics::StateVectorObserver);
%fp_shared_ptr(FullPhysics::SubStateVectorObserver);
namespace FullPhysics {
class StateVector;
}

%fp_shared_ptr(FullPhysics::Observer<FullPhysics::StateVector>);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::StateVector>);

// Do this so we can derive from this and have it able to be used by the C++ code
// Defined here since rename does not like being inside of a namespace
%feature("director") FullPhysics::Observer<FullPhysics::StateVector>;
%rename(ObserverStateVector) FullPhysics::Observer<FullPhysics::StateVector>;

namespace FullPhysics {
class FullPhysics::Observer<FullPhysics::StateVector> {
public:
  virtual ~Observer();
  virtual void notify_add(StateVector& Obs);
  virtual void notify_remove(StateVector& Obs);
  virtual void notify_update(const StateVector& Obs);
};
}

namespace FullPhysics {
%template(ObservableStateVector) FullPhysics::Observable<FullPhysics::StateVector>;

%nodefaultctor StateVectorObserver;

class StateVectorObserver : public Observer<StateVector> {
public:
  virtual ~StateVectorObserver();
  std::string print_to_string() const;
  virtual void mark_used(const StateVector& Sv, 
			 blitz::Array<bool, 1>& Used) const;
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;
};

class StateVector: public Observable<StateVector> {
public:
  virtual ~StateVector();
  std::string print_to_string() const;
  StateVector();
  virtual void add_observer(Observer<StateVector>& Obs);
  virtual void remove_observer(Observer<StateVector>& Obs);
  %python_attribute(state, blitz::Array<double, 1>);
  %python_attribute(state_with_derivative, ArrayAd<double, 1>);
  %extend {
    std::vector<std::string> _state_vector_name() const 
    {
      blitz::Array<std::string, 1> sn($self->state_vector_name());
      std::vector<std::string> res;
      for(int i = 0; i < sn.extent(blitz::firstDim); ++i)
	res.push_back(sn(i));
      return res;
    }
  }
%pythoncode {
@property
def state_vector_name(self):
    return self._state_vector_name()
}
  %python_attribute(state_covariance, blitz::Array<double, 2>);
  void update_state(const blitz::Array<double, 1>& X);
  void update_state(const blitz::Array<double, 1>& X, 
		    const blitz::Array<double, 2>& Cov);
  %python_attribute(used_flag, blitz::Array<bool, 1>);
  %python_attribute_with_set(observer_claimed_size, int);
};

class SubStateVectorObserver : public StateVectorObserver {
public:
  virtual ~SubStateVectorObserver();
  virtual void notify_update(const StateVector& Sv);
  virtual void mark_used(const StateVector& Sv, 
			 blitz::Array<bool, 1>& Used) const;
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;
  %python_attribute(state_vector_start_index, int);
  %python_attribute(sub_vector_size, int);
  virtual void update_sub_state(
    const ArrayAd<double, 1>& Sv_sub,
    const blitz::Array<double, 2>& Cov_sub) = 0;
  virtual void mark_used_sub(blitz::Array<bool, 1>& Used) const;
  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name)
    const;
  virtual void print(std::ostream& Os) const;
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
protected:
  SubStateVectorObserver(int Plen);
  SubStateVectorObserver();
  void state_vector_observer_initialize(int Plen);
  blitz::Array<double, 1> sv_full;
  blitz::Array<double, 2> sv_cov_full;
  blitz::Array<double, 1> sv_sub;
  blitz::Array<double, 2> sv_cov_sub;
};

}
