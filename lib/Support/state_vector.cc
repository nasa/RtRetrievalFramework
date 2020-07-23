#include "state_vector.h"
#include "spectrum_effect.h"
#include "fp_exception.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
#include "rt_atmosphere.h"
#include "instrument.h"
#include "stokes_coefficient.h"

void state_vector_add_observer_instrument(StateVector& Sv, Instrument& inst)
{
  Sv.add_observer(inst);
}

void state_vector_add_observer_atm(StateVector& Sv, RtAtmosphere& atm)
{
  Sv.add_observer(atm);
}

void state_vector_add_observer_scoeff(StateVector& Sv, StokesCoefficient& Scoef)
{
  Sv.add_observer(Scoef);
}

void state_vector_add_observer_spec_effect(StateVector& Sv, std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& se)
{
  BOOST_FOREACH(std::vector<boost::shared_ptr<SpectrumEffect> >& i, se) {
    BOOST_FOREACH(boost::shared_ptr<SpectrumEffect>& j, i)
      Sv.add_observer(*j);
  }
}

void state_vector_print_names(StateVector& Sv)
{
  BOOST_FOREACH(const std::string& s, Sv.state_vector_name())
    std::cout << s << "\n";
}
typedef void (StateVector::*us1)(const blitz::Array<double, 1>&);
typedef void (StateVector::*us2)(const blitz::Array<double, 1>&, const blitz::Array<double, 2>&);
REGISTER_LUA_CLASS(StateVector)
.def(luabind::constructor<>())
.def("add_observer", &state_vector_add_observer_instrument)
.def("add_observer", &state_vector_add_observer_atm)
.def("add_observer", &state_vector_add_observer_scoeff)
.def("add_observer", &state_vector_add_observer_spec_effect)
.def("update_state", ((us1) &StateVector::update_state))
.def("update_state", ((us2) &StateVector::update_state))
.def("print_names", &state_vector_print_names)
.def("state", &StateVector::state)
.def("state_covariance", &StateVector::state_covariance)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Update the state vector. This version sets the covariance to a
/// dummy identity matrix that is the same size as state().
//-----------------------------------------------------------------------

void StateVector::update_state(const blitz::Array<double, 1>& X)
{
  x_.resize(X.rows(), X.rows());
  x_.value() = X;
  x_.jacobian() = 0;
  for(int i = 0; i < x_.rows(); ++i)
    x_.jacobian()(i,i) = 1;
  cov_.resize(x_.rows(), x_.rows());
  cov_ = 0;
  for(int i = 0; i < x_.rows(); ++i)
    cov_(i, i) = 1;
  notify_update_do(*this);
}

//-----------------------------------------------------------------------
/// Update the state vector and covariance.
//-----------------------------------------------------------------------

void StateVector::update_state(const blitz::Array<double, 1>& X, 
			       const blitz::Array<double, 2>& Cov)
{
  if(X.rows() != Cov.rows() ||
     X.rows() != Cov.cols())
    throw Exception("X and Cov need to be the same size when updating the StateVector");
  x_.resize(X.rows(), X.rows());
  x_.value() = X;
  x_.jacobian() = 0;
  for(int i = 0; i < x_.rows(); ++i)
    x_.jacobian()(i,i) = 1;
  cov_.reference(Cov.copy());
  notify_update_do(*this);
}


//-----------------------------------------------------------------------
/// Return a Array of boolean values. The value (i) is true if the state
/// vector element X(i) is being used. This can be used to determine
/// parameters that are being ignored, e.g. the number of active
/// levels in an Aerosol is less that the size of the state vector for it.
//-----------------------------------------------------------------------

blitz::Array<bool, 1> StateVector::used_flag() const
{
  blitz::Array<bool, 1> res(state().rows());
  res = false;
  BOOST_FOREACH(const boost::weak_ptr<Observer<StateVector> >& t, olist) {
    boost::shared_ptr<StateVectorObserver> t2 = 
      boost::dynamic_pointer_cast<StateVectorObserver>(t.lock());
    if(t2)
      t2->mark_used(*this, res);
  }
  return res;
}

void SubStateVectorObserver::notify_update(const StateVector& Sv)
{
  if(Sv.state().rows() < pstart + plen) {
    // blitz::Array<std::string, 1> svname =Sv.state_vector_name();
    // for(int i = 0; i < svname.rows(); ++i) {
    //   if(i < Sv.state().rows())
    // 	std::cerr << Sv.state()(i) << " : " << svname(i) << "\n";
    //   else
    // 	std::cerr <<  "Out_of_range : " << svname(i) << "\n";
    // }
    std::stringstream err_msg;
    err_msg << "StateVector is of size: "
	    << Sv.state().rows()
	    << " not the expected size: "
	    << (pstart + plen);
    throw Exception(err_msg.str());
  }
  if(pstart < 0)
    throw Exception("pstart < 0");
  sv_full.reference(Sv.state_with_derivative());
  sv_cov_full.reference(Sv.state_covariance());
  if(plen > 0) {
    blitz::Range rsub(pstart, pstart + plen - 1);
    sv_sub.reference(sv_full(rsub));
    sv_cov_sub.reference(Array<double, 2>(sv_cov_full(rsub, rsub)));
  }
  // So that anything registered as a SubStateVector observer
  // gets notified even if it doesn't contain any SV elements
  update_sub_state(sv_sub, sv_cov_sub);
}

void SubStateVectorObserver::mark_used(const StateVector& Sv, 
				blitz::Array<bool, 1>& Used) const
{
  if(Used.rows() < pstart + plen)
    throw Exception("StateVector not the expected size");
  if(pstart < 0)
    throw Exception("pstart < 0");
  if(plen > 0) {
    Array<bool, 1> used_sub(Used(blitz::Range(pstart, pstart + plen - 1)));
    mark_used_sub(used_sub);
  }
}

void SubStateVectorObserver::state_vector_name(const StateVector& Sv, 
		       blitz::Array<std::string, 1>& Sv_name) const
{
  if(Sv_name.rows() < pstart + plen)
    throw Exception("StateVector not the expected size");
  if(pstart < 0)
    throw Exception("pstart < 0");
  if(plen > 0) {
    Array<std::string, 1> sv_name_sub(Sv_name(blitz::Range(pstart, 
							   pstart + plen - 1)));
    state_vector_name_sub(sv_name_sub);
  }
}

//-----------------------------------------------------------------------
/// Return name of each state vector element.
//-----------------------------------------------------------------------

Array<std::string, 1> StateVector::state_vector_name() const
{
  Array<std::string, 1> res(std::max(state().rows(), observer_claimed_size()));
  for(int i = 0; i < res.rows(); ++i)
    res(i) = "State vector " + boost::lexical_cast<std::string>(i + 1);
  BOOST_FOREACH(const boost::weak_ptr<Observer<StateVector> >& t, olist) {
    boost::shared_ptr<StateVectorObserver> t2 = 
      boost::dynamic_pointer_cast<StateVectorObserver>(t.lock());
    if(t2)
      t2->state_vector_name(*this, res);
  }
  return res;
}

void StateVector::print(std::ostream& Os) const
{
  const static int sv_num_width = 17;
  Array<std::string, 1> svname = state_vector_name();
  for(int i = 0; i < std::max(svname.rows(), state().rows()); ++i) {
    Os << std::setprecision(sv_num_width-7)
       << std::setw(sv_num_width);

    if(i < state().rows())
      Os << state()(i);
    else
      Os << "Past end state vector";

    Os << "  ";
    if(i < svname.rows())
      Os << svname(i);
    else
      Os << "Unlabeled row " << i;
    Os << std::endl;

  }
}

//-----------------------------------------------------------------------
/// Take the given number of state vector parameters. We determine
/// where the starting point to use is when we attach to the state
/// vector. 
///
/// Note that it is perfectly legal for Plen to be 0, that just means
/// we don't have any parameters. This is a useful edge case that we
/// support. 
//-----------------------------------------------------------------------

void SubStateVectorObserver::state_vector_observer_initialize(int Plen)
{
  range_min_check(Plen, 0);
  plen = Plen;
}
