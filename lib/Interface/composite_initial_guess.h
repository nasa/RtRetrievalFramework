#ifndef COMPOSITE_INITIAL_GUESS_H
#define COMPOSITE_INITIAL_GUESS_H
#include "initial_guess.h"
#include <vector>
#include <list>
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
class CompositeInitialGuess;
/****************************************************************//**
  Class that builds a portion of the state vector.

  We use a std::vector here rather than a blitz::Array just because it
  is easier to add something to the end of std::vector than
  blitz::Array. CompositeInitialGuess converts this to a blitz::Array
  before finishing the construction of the initial guess.
*******************************************************************/

class InitialGuessBuilder : public Printable<InitialGuessBuilder> {
public:
  virtual ~InitialGuessBuilder() {}

//-----------------------------------------------------------------------
/// Called when we get attached to a CompositeInitialGuess. The default 
/// is to do nothing, but derived classes can override this if desired.
//-----------------------------------------------------------------------

  virtual void attach_notify(CompositeInitialGuess& Comp_ig)
  { }

//-----------------------------------------------------------------------
/// Number of elements we will be adding to the state vector. 0 is a legal
/// value, if we are changing elements but not adding any.
//-----------------------------------------------------------------------

  virtual int number_element() const = 0;

//-----------------------------------------------------------------------
/// Called when we need this class to do its part in setting up the
/// initial state vector. 
///
/// \param v State vector that should be updated in place.
/// \param index Since we are often adding to the end of the state
///     vector, index is passed in. This is the sum of the
///     number_elements() of all the 
///     InitialGuessBuilder that appear before this object in the list.
//-----------------------------------------------------------------------  

  virtual void build_initial_value(blitz::Array<double, 1>& v, int index) 
    const = 0;

//-----------------------------------------------------------------------
/// Called when we need this class to do its part in setting up the
/// apriori state vector. 
///
/// \param v State vector that should be updated in place.
/// \param index Since we are often adding to the end of the state
///     vector, index is passed in. This is the sum of the
///     number_elements() of all the 
///     InitialGuessBuilder that appear before this object in the list.
//-----------------------------------------------------------------------  

  virtual void build_apriori(blitz::Array<double, 1>& v, int index) const = 0;

//-----------------------------------------------------------------------
/// Called when we need this class to do its part in setting up the
/// covariance matrix for the a priori state
/// vector. 
///
/// \param m State vector that should be updated in place.
/// \param index Since we are often adding to the end of the state
///     vector, index is passed in. This is the sum of the
///     number_elements() of all the 
///     InitialGuessBuilder that appear before this object in the list.
//-----------------------------------------------------------------------  

  virtual void build_apriori_covariance(blitz::Array<double, 2>& m, 
					int index) const = 0;
  virtual void print(std::ostream& Os) const {Os << "InitialGuessBuilder";}
};

/****************************************************************//**
  A common way to create an initial guess is to have other classes
  responsible for portions of the state vector (e.g., an Atmosphere
  class creates the portion of the initial guess that handles the
  description of the atmosphere layers). This class implements this
  division. 

  This is an example of the "Builder" design pattern. This class is
  what is commonly called the "Director", and the InitialGuessBuilder
  classes are the "Builder" classes.

  Note that the InitialGuessBuilder objects are called in the order
  they are added to the CompositeInitialGuess object. A common
  InitialGuessBuilder adds additional values to the end state vector, so
  the order is important.

  A CompositeInitialGuess is also a InitialGuessBuilder, so you can
  use this to nest initial guess builders.
*******************************************************************/

class CompositeInitialGuess : public InitialGuess,
			      public InitialGuessBuilder {
public:
  virtual ~CompositeInitialGuess() {}
  virtual blitz::Array<double, 1> initial_guess() const;
  virtual blitz::Array<double, 1> apriori() const;
  virtual blitz::Array<double, 2> apriori_covariance() const;
  virtual void print(std::ostream& Os) const;
  virtual int number_element() const;

  virtual void build_initial_value(blitz::Array<double, 1>& v, int index) 
    const;
  virtual void build_apriori(blitz::Array<double, 1>& v, int index) const;
  virtual void build_apriori_covariance(blitz::Array<double, 2>& m, 
					int index) const;

//-----------------------------------------------------------------------
/// Add a builder to the build list.
//-----------------------------------------------------------------------

  void add_builder(const boost::shared_ptr<InitialGuessBuilder>& B)
  {
    blist.push_back(B);
    B->attach_notify(*this);
  }

//-----------------------------------------------------------------------
/// Remove a builder to the build list.
//-----------------------------------------------------------------------

  void remove_builder(const boost::shared_ptr<InitialGuessBuilder>& B)
  {
    blist.remove(B);
  }
private:
  std::list<boost::shared_ptr<InitialGuessBuilder> > blist;
};
}
#endif
