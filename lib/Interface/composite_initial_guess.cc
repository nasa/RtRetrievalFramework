#include "composite_initial_guess.h"
#include <boost/foreach.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua has trouble with the multiple inheritance of
// CompositeInitialGuess, so add a add_builder that has the right type
// without needed to convert to InitialGuessBuilder
void composite_initial_guess_add_builder(CompositeInitialGuess& Ig, 
 const boost::shared_ptr<InitialGuess>& V)
{
  Ig.add_builder(boost::dynamic_pointer_cast<CompositeInitialGuess>(V));
}
void composite_initial_guess_add_builder2(CompositeInitialGuess& Ig, 
 const boost::shared_ptr<CompositeInitialGuess>& V)
{
  Ig.add_builder(boost::dynamic_pointer_cast<CompositeInitialGuess>(V));
}

REGISTER_LUA_CLASS(InitialGuessBuilder)
REGISTER_LUA_END()
REGISTER_LUA_DERIVED_CLASS(CompositeInitialGuess, InitialGuess)
.def(luabind::constructor<>())
.def("add_builder", &CompositeInitialGuess::add_builder)
.def("add_builder", &composite_initial_guess_add_builder) 
.def("add_builder", &composite_initial_guess_add_builder2) 
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Return the initial state vector to use.
//-----------------------------------------------------------------------

blitz::Array<double, 1> CompositeInitialGuess::initial_guess() const
{
  blitz::Array<double, 1> res(number_element());
  int index = 0;
  BOOST_FOREACH(const boost::shared_ptr<InitialGuessBuilder>& v, blist) {
    v->build_initial_value(res, index);
    index += v->number_element();
  }
  return res;
}

//-----------------------------------------------------------------------
/// Number of elements that will be in the state vector.
//-----------------------------------------------------------------------

int CompositeInitialGuess::number_element() const
{
  int res = 0;
  BOOST_FOREACH(const boost::shared_ptr<InitialGuessBuilder>& v, blist)
    res += v->number_element();
  return res;
}

//-----------------------------------------------------------------------
/// Return the apriori state vector to use.
//-----------------------------------------------------------------------

blitz::Array<double, 1> CompositeInitialGuess::apriori() const
{
  blitz::Array<double, 1> res(number_element());
  int index = 0;
  BOOST_FOREACH(const boost::shared_ptr<InitialGuessBuilder>& v, blist) {
    v->build_apriori(res, index);
    index += v->number_element();
  }
  return res;
}

//-----------------------------------------------------------------------
/// Return the apriori state vector covariance to use.
//-----------------------------------------------------------------------

blitz::Array<double, 2> CompositeInitialGuess::apriori_covariance() 
const
{
  int n = number_element();
  // Create array in Fortran Column major order, but still use a 0 base.
  blitz::Array<double, 2> res(n,n,blitz::ColumnMajorArray<2>());
  res = 0.0;
  int index = 0;
  BOOST_FOREACH(const boost::shared_ptr<InitialGuessBuilder>& v, blist) {
    v->build_apriori_covariance(res, index);
    index += v->number_element();
  }
  return res;
}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void CompositeInitialGuess::print(std::ostream& Os) const
{
  Os << "Composite Initial Guess";
}

void CompositeInitialGuess::build_initial_value
(blitz::Array<double, 1>& v, int index) const
{
  if(number_element() == 0)
    return;
  v(Range(index, index + number_element() - 1)) = initial_guess();
}

void CompositeInitialGuess::build_apriori
(blitz::Array<double, 1>& v, int index) const
{
  if(number_element() == 0)
    return;
  v(Range(index, index + number_element() - 1)) = apriori();
}

void CompositeInitialGuess::build_apriori_covariance
(blitz::Array<double, 2>& m, int index) const
{
  if(number_element() == 0)
    return;
  Range r(index, index + number_element() - 1);
  m(r, r) = apriori_covariance();
}



