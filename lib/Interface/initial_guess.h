#ifndef INITIAL_GUESS_H
#define INITIAL_GUESS_H
#include "printable.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This gets the initial guess and the apriori state vector
  values. Often but not always the initial guess to use is the same as
  the a priori.
*******************************************************************/
class InitialGuess : public Printable<InitialGuess> {
public:
  virtual ~InitialGuess() {}

//-----------------------------------------------------------------------
/// Return the initial state vector to use.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> initial_guess() const = 0;

//-----------------------------------------------------------------------
/// Return the apriori state vector to use.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> apriori() const = 0;

//-----------------------------------------------------------------------
/// Return the apriori state vector covariance to use.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> apriori_covariance() const = 0;

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const {Os << "InitialGuess";}
};
}
#endif
