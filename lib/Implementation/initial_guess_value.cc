#include "initial_guess_value.h"
#include "ostream_pad.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
typedef const blitz::Array<double, 1>& (InitialGuessValue::*a1)(void) const;
typedef void (InitialGuessValue::*a2)(const blitz::Array<double, 1>&);
typedef const blitz::Array<double, 2>& (InitialGuessValue::*b1)(void) const;
typedef void (InitialGuessValue::*b2)(const blitz::Array<double, 2>&);
double initial_guess_value_apriori_double_get(const InitialGuessValue& ig)
{ return ig.apriori()(0); }
void initial_guess_value_apriori_double_set(InitialGuessValue& ig, 
					      double val)
{ Array<double, 1> v(1);
  v(0) = val;
  ig.apriori(v);
}
double initial_guess_value_initial_guess_double_get(const InitialGuessValue& ig)
{ return ig.initial_guess()(0); }
void initial_guess_value_initial_guess_double_set(InitialGuessValue& ig, 
					      double val)
{ Array<double, 1> v(1);
  v(0) = val;
  ig.initial_guess(v);
}
double initial_guess_value_apriori_cov_double_get(const InitialGuessValue& ig)
{ return ig.apriori_covariance()(0, 0); }
void initial_guess_value_apriori_cov_double_set(InitialGuessValue& ig, 
					      double val)
{ Array<double, 2> v(1,1);
  v(0,0) = val;
  ig.apriori_covariance(v);
}

REGISTER_LUA_DERIVED_CLASS(InitialGuessValue, InitialGuessBuilder)
.def(luabind::constructor<>())
.property("apriori", 
	  ((a1) &InitialGuessValue::apriori),
	  ((a2) &InitialGuessValue::apriori))
.def("apriori_subset", &InitialGuessValue::apriori_subset)
.def("initial_guess_subset", &InitialGuessValue::initial_guess_subset)
.def("apriori_covariance_subset", 
     &InitialGuessValue::apriori_covariance_subset)
.property("apriori_double", 
	  &initial_guess_value_apriori_double_get,
	  &initial_guess_value_apriori_double_set)
.property("initial_guess", 
	  ((a1) &InitialGuessValue::initial_guess),
	  ((a2) &InitialGuessValue::initial_guess))
.property("initial_guess_double", 
	  &initial_guess_value_initial_guess_double_get,
	  &initial_guess_value_initial_guess_double_set)
.property("apriori_covariance", 
	  ((b1) &InitialGuessValue::apriori_covariance),
	  ((b2) &InitialGuessValue::apriori_covariance))
.property("apriori_covariance_double", 
	  &initial_guess_value_apriori_cov_double_get,
	  &initial_guess_value_apriori_cov_double_set)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Subset a value to include only those elements where Flag is true.
//-----------------------------------------------------------------------

void InitialGuessValue::apriori_subset
(const blitz::Array<bool, 1>& Flag, 
 const blitz::Array<double, 1>& V)
{
  if(Flag.rows() != V.rows()) {
    std::stringstream err_msg;
    err_msg << "Flag and apriori need to be same size. "
	    << Flag.rows() << " != " << V.rows();
    throw Exception(err_msg.str());
  }
  apriori_.resize(count(Flag));
  int ind = 0;
  for(int i = 0; i < V.rows(); ++i)
    if(Flag(i))
      apriori_(ind++) = V(i);
}

//-----------------------------------------------------------------------
/// Subset a value to include only those elements where Flag is true.
//-----------------------------------------------------------------------

void InitialGuessValue::initial_guess_subset
(const blitz::Array<bool, 1>& Flag, 
 const blitz::Array<double, 1>& V)
{
  if(Flag.rows() != V.rows()) {
    std::stringstream err_msg;
    err_msg << "Flag and initial guess need to be same size. "
            << Flag.rows() << " != " << V.rows();
    throw Exception(err_msg.str());
  }
  initial_guess_.resize(count(Flag));
  int ind = 0;
  for(int i = 0; i < V.rows(); ++i)
    if(Flag(i))
      initial_guess_(ind++) = V(i);
}

//-----------------------------------------------------------------------
/// Subset a value to include only those elements where Flag is true.
//-----------------------------------------------------------------------

void InitialGuessValue::apriori_covariance_subset
(const blitz::Array<bool, 1>& Flag, 
 const blitz::Array<double, 2>& V)
{
  if(Flag.rows() != V.rows() ||
     Flag.rows() != V.cols()) {
    std::stringstream err_msg;
    err_msg << "Flag and covariance need to be same size. "
	    << Flag.rows() << " != " << "[ " << V.rows() << ", " << V.cols() << " ]";
    throw Exception(err_msg.str());
  }

  apriori_covariance_.resize(count(Flag), count(Flag));
  int ind1 = 0;
  for(int i = 0; i < V.rows(); ++i)
    if(Flag(i)) {
      int ind2 = 0;
      for(int j = 0; j < V.cols(); ++j)
	if(Flag(j))
	  apriori_covariance_(ind1, ind2++) = V(i, j);
      ++ind1;
    }
}

void InitialGuessValue::build_initial_value(blitz::Array<double, 1>& v, 
					    int index) const
{
  if(number_element() == 0)
    return;
  if(initial_guess().rows() != number_element()) {
    Exception e;
    e << "initial_guess() in InitialGuessValue has " << initial_guess().rows()
      << " elements, but apriori has " << number_element();
    throw e;
  }
  v(Range(index, index + number_element() - 1)) = initial_guess();
}

void InitialGuessValue::build_apriori(blitz::Array<double, 1>& v, 
				      int index) const
{
  if(number_element() == 0)
    return;
  v(Range(index, index + number_element() - 1)) = apriori();
}

void InitialGuessValue::build_apriori_covariance(blitz::Array<double, 2>& m, 
						 int index) const
{
  if(number_element() == 0)
    return;
  if(apriori_covariance().rows() != number_element() ||
     apriori_covariance().cols() != number_element()) {
    Exception e;
    e << "apriori_covariance() in InitialGuessValue has " 
      << apriori_covariance().rows() << " x " << apriori_covariance().cols()
      << " elements, but apriori has " << number_element();
    throw e;
  }
  Range r(index, index + number_element() - 1);
  m(r, r) = apriori_covariance();
}

void InitialGuessValue::print(std::ostream& Os) const 
{
  Os << "InitialGuessValue:\n";
  OstreamPad opad(Os, "  ");
  Os << "Apriori:\n" << apriori() << "\n"
     << "Initial guess:\n" << initial_guess() << "\n"
     << "Covariance:\n" << apriori_covariance() << "\n";
}
