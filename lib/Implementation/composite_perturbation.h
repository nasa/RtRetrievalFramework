#ifndef COMPOSITE_PERTURBATION_H
#define COMPOSITE_PERTURBATION_H
#include "perturbation.h"
#include <vector>
#include <list>
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
class CompositePerturbation;
/****************************************************************//**
  Class that builds a perturbation to use for a finite difference
  Jacobian. 

  We use a std::vector here rather than a blitz::Array just because it
  is easier to add something to the end of std::vector than
  blitz::Array. CompositePerturbation converts this to a blitz::Array
  before finishing the construction of the initial guess.
*******************************************************************/

class PerturbationBuilder {
public:
  virtual ~PerturbationBuilder() {}

//-----------------------------------------------------------------------
/// Called when we get attached to a CompositePerturbation. The default 
/// is to do nothing, but derived classes can override this if desired.
//-----------------------------------------------------------------------

  virtual void attach_notify(CompositePerturbation& Comp_ig)
  { }

//-----------------------------------------------------------------------
/// Number of elements we will be adding to the perturbation. 0 is a
/// legal value, if we are changing elements but not adding any.
//-----------------------------------------------------------------------

  virtual int number_element() const = 0;

//-----------------------------------------------------------------------
/// Called when we need this class to do its part in setting up the
/// perturbation array.
///
/// \param v Perturbation vector that should be updated in place.
/// \param index Since we are often adding to the end of the state
///     vector, index is passed in. This is the sum of the
///     number_elements() of all the 
///     PerturbationBuilder that appear before this object in the list.
//-----------------------------------------------------------------------  

  virtual void build_perturbation(blitz::Array<double, 1>& v, int index) const = 0;
};

/****************************************************************//**
  A common way to create a perturbation is to have other classes
  responsible for portions of the state vector (e.g., an Atmosphere
  class creates the portion of the initial guess that handles the
  description of the atmosphere layers). This class implements this
  division. 

  This is an example of the "Builder" design pattern. This class is
  what is commonly called the "Director", and the PerturbationBuilder
  classes are the "Builder" classes.

  Note that the PerturbationBuilder objects are called in the order
  they are added to the CompositePerturbation object. A common
  PerturbationBuilder adds additional values to the end state vector, so
  the order is important.
*******************************************************************/

class CompositePerturbation : public Perturbation {
public:
  virtual ~CompositePerturbation() {}
  virtual blitz::Array<double, 1> perturbation() const;
  virtual void print(std::ostream& Os) const;
  int number_element() const;
//-----------------------------------------------------------------------
/// Add a builder to the build list.
//-----------------------------------------------------------------------

  void add_builder(const boost::shared_ptr<PerturbationBuilder>& B)
  {
    blist.push_back(B);
    B->attach_notify(*this);
  }

//-----------------------------------------------------------------------
/// Remove a builder to the build list.
//-----------------------------------------------------------------------

  void remove_builder(const boost::shared_ptr<PerturbationBuilder>& B)
  {
    blist.remove(B);
  }
private:
  std::list<boost::shared_ptr<PerturbationBuilder> > blist;
};
}
#endif
