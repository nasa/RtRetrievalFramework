#ifndef STATE_VECTOR_H
#define STATE_VECTOR_H
#include "printable.h"
#include "observer.h"
#include "array_ad.h"
#include "fp_exception.h"
#include <blitz/array.h>

namespace FullPhysics {
class StateVector;
/****************************************************************//**
  This is an observer of a StateVector. If attached to a StateVector,
  this class gets notified when a state vector is updated.

  It is completely unspecified what an observer does with this
  information, but commonly the class will update its internal state
  based on the state vector update.
*******************************************************************/
class StateVectorObserver : public Printable<StateVectorObserver>,
  public Observer<StateVector> {
public:
  virtual ~StateVectorObserver() {}

//-----------------------------------------------------------------------
/// Mark elements that we are actively using (i.e., that aren't
/// ignored). You only need to mark the ones that are used as true,
/// everything is already initialized as false. Default is to do nothing.
//-----------------------------------------------------------------------

  virtual void mark_used(const StateVector& Sv, 
			 blitz::Array<bool, 1>& Used) const {}
			 
//-----------------------------------------------------------------------
/// Update any portion of the list of the state vector names that
/// apply to this object. Default is to do nothing.
//-----------------------------------------------------------------------

  virtual void state_vector_name(const StateVector& Sv, 
				 blitz::Array<std::string, 1>& Sv_name) const {}
			 
  virtual void print(std::ostream& Os) const { Os << "StateVectorObserver";}
};

/****************************************************************//**
  This handles informing a set of interested objects when the state
  vector has updated. Those objects then update their internal state
  to account for the new state vector.
*******************************************************************/
class StateVector : public Printable<StateVector>, 
  public Observable<StateVector> {
public:
  StateVector()
    : pstart(0)
  {
  }
  virtual ~StateVector() { }

  virtual void add_observer(Observer<StateVector>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<StateVector>& Obs) 
  { remove_observer_do(Obs, *this);}
  virtual void print(std::ostream& Os) const;

//-----------------------------------------------------------------------
/// Current state vector.
//-----------------------------------------------------------------------

  const blitz::Array<double, 1>& state() const { return x_.value(); }

//-----------------------------------------------------------------------
/// Return the state vector as state() does, but also make each value
/// a AutoDerivative. The derivative is with respect to the state
/// vector, i.e., we treat the state vector as the independent
/// variables. This means the first value has a gradient all 0's
/// except for 1 in the first index, the second value all zeros except
/// for 1 in the second index, etc.
//-----------------------------------------------------------------------

  const ArrayAd<double, 1>& state_with_derivative() const { return x_;}

  blitz::Array<std::string, 1> state_vector_name() const;

//-----------------------------------------------------------------------
/// Current covariance of the state vector.
//-----------------------------------------------------------------------

  const blitz::Array<double, 2>& state_covariance() const {return cov_;}

  void update_state(const blitz::Array<double, 1>& X);
  void update_state(const blitz::Array<double, 1>& X, 
		    const blitz::Array<double, 2>& Cov);
  blitz::Array<bool, 1> used_flag() const;

//-----------------------------------------------------------------------
/// Total "claimed" size of the state vector. For observers that
/// register an interest in a portion of the state vector, we add all
/// of the portions of interest. Note that an actual state vector
/// isn't constrained to this size, it might be larger (with
/// presumably portions ignored), or if the observers handle it
/// correctly it could be smaller.
//-----------------------------------------------------------------------

  int observer_claimed_size() const {return pstart;}

//-----------------------------------------------------------------------
/// Update claimed size of state vector.
//-----------------------------------------------------------------------

  void observer_claimed_size(int Pstart) { pstart = Pstart; }

private:
  ArrayAd<double, 1> x_;
  blitz::Array<double, 2> cov_;
  // Helper value that says what portion of state vector has been
  // claimed by observers. We don't do anything with this value in
  // this class, except make it available to StateVectorObservers when 
  // they are attached.
  int pstart;			
};

/****************************************************************//**
  A common StateVectorObserver just "owns" a subset of the
  StateVector. This class gives the common behavior for this case.
*******************************************************************/
class SubStateVectorObserver : virtual public StateVectorObserver {
public:
  virtual ~SubStateVectorObserver() {}
  virtual void notify_update(const StateVector& Sv);
  virtual void mark_used(const StateVector& Sv, 
			 blitz::Array<bool, 1>& Used) const;
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;

//-----------------------------------------------------------------------
/// Starting index of state vector used by this object.
//-----------------------------------------------------------------------

  int state_vector_start_index() const {return pstart; }

//-----------------------------------------------------------------------
/// Length of the sub set of the state vector used by this object.
//-----------------------------------------------------------------------

  int sub_vector_size() const {return plen; }

//-----------------------------------------------------------------------
/// Called by update_state with the subset of the state vector used by
/// this class. 
//-----------------------------------------------------------------------

  virtual void update_sub_state(
    const ArrayAd<double, 1>& Sv_sub,
    const blitz::Array<double, 2>& Cov_sub) = 0;

//-----------------------------------------------------------------------
/// Called by mark_used with the subset of the state vector used by
/// this class. The default marks everything as used, but derived
/// classes can override this.
//-----------------------------------------------------------------------

  virtual void mark_used_sub(blitz::Array<bool, 1>& Used) const
  { Used = true; }

//-----------------------------------------------------------------------
/// Called by state_vector_name with the subset of the Sv_name used by
/// this class. The default function doesn't change anything, but
/// derived classes can ovveride this.
//-----------------------------------------------------------------------

  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) 
    const {}
  virtual void print(std::ostream& Os) const {Os << "SubStatVectorObserver";}
  virtual void notify_add(StateVector& Sv)
  { 
    if(pstart != -1)
      throw Exception("A SubStateVectorObserver can only be attached to one state vector");
    pstart = Sv.observer_claimed_size();
    Sv.observer_claimed_size(pstart + plen);
  }

  virtual void notify_remove(StateVector& Sv)
  {
    pstart = -1;
  }
protected:
//-----------------------------------------------------------------------
/// Take the given number of state vector parameters. We determine
/// where the starting point to use is when we attach to the state
/// vector. 
///
/// Note that it is perfectly legal for Plen to be 0, that just means
/// we don't have any parameters. This is a useful edge case that we
/// support. 
//-----------------------------------------------------------------------

  SubStateVectorObserver(int Plen)
  : pstart(-1)
  {
    state_vector_observer_initialize(Plen);
  }

//-----------------------------------------------------------------------
/// Default constructor. Derived classes should call initialize before
/// finishing there constructor.
//-----------------------------------------------------------------------

  SubStateVectorObserver() : pstart(-1), plen(0) {}

  void state_vector_observer_initialize(int Plen);

//-----------------------------------------------------------------------
/// The last full state vector we have been updated with, saved for
/// reference by derived class
//-----------------------------------------------------------------------
  ArrayAd<double, 1> sv_full;

//-----------------------------------------------------------------------
/// The last full covariance matrix we have been with, saved for
/// reference by derived class.
//-----------------------------------------------------------------------
  blitz::Array<double, 2> sv_cov_full;

//-----------------------------------------------------------------------
/// The subset of sv_full that is "owned" by this class, what was
/// passed through update_sub_state. Saved for reference by derived
/// class. 
//-----------------------------------------------------------------------
  ArrayAd<double, 1> sv_sub;

//-----------------------------------------------------------------------
/// The subset of cov_full that is "owned" by this class, what was
/// passed through update_sub_state. Saved for reference by derived
/// class. 
//-----------------------------------------------------------------------
  blitz::Array<double, 2> sv_cov_sub;
private:
  int pstart;
  int plen;
};
}
#endif
