#include "composite_perturbation.h"
#include <boost/foreach.hpp>

using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Number of elements that will be in the state vector.
//-----------------------------------------------------------------------

int CompositePerturbation::number_element() const
{
  int res = 0;
  BOOST_FOREACH(const boost::shared_ptr<PerturbationBuilder>& v, blist)
    res += v->number_element();
  return res;
}

//-----------------------------------------------------------------------
/// Return the perturbation vector to use.
//-----------------------------------------------------------------------

blitz::Array<double, 1> CompositePerturbation::perturbation() const
{
  blitz::Array<double, 1> res(number_element());
  int index = 0;
  BOOST_FOREACH(const boost::shared_ptr<PerturbationBuilder>& v, blist) {
    v->build_perturbation(res, index);
    index += v->number_element();
  }
  return res;
}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void CompositePerturbation::print(std::ostream& Os) const
{
  Os << "Composite Perturbation";
}

