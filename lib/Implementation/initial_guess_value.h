#ifndef INITIAL_GUESS_VALUE_H
#define INITIAL_GUESS_VALUE_H
#include "composite_initial_guess.h"

namespace FullPhysics {
/****************************************************************//**
  This is a simple implementation of InitialGuessBuilder that just
  has variables used to give the apriori, initial guess, and
  covariance. Since the initial guess of often just the apriori, the
  initial guess is allowed to be left empty and we just use the apriori
  in that case.
*******************************************************************/

class InitialGuessValue : public InitialGuessBuilder {
public:
  virtual ~InitialGuessValue() {}
  virtual int number_element() const { return apriori_.rows(); }
  virtual void build_initial_value(blitz::Array<double, 1>& v, int index) 
    const;
  virtual void build_apriori(blitz::Array<double, 1>& v, int index) const;
  virtual void build_apriori_covariance(blitz::Array<double, 2>& m, 
					int index) const;

//-----------------------------------------------------------------------
/// Apriori value
//-----------------------------------------------------------------------

  const blitz::Array<double, 1>& apriori() const {return apriori_;}

//-----------------------------------------------------------------------
/// Set apriori value
//-----------------------------------------------------------------------

  void apriori(const blitz::Array<double, 1>& v) {apriori_.reference(v.copy());}

  void apriori_subset(const blitz::Array<bool, 1>& Flag, 
		      const blitz::Array<double, 1>& V);
  void initial_guess_subset(const blitz::Array<bool, 1>& Flag, 
                            const blitz::Array<double, 1>& V);
  void apriori_covariance_subset(const blitz::Array<bool, 1>& Flag, 
				 const blitz::Array<double, 2>& V);

//-----------------------------------------------------------------------
/// Apriori covariance value
//-----------------------------------------------------------------------

  const blitz::Array<double, 2>& apriori_covariance() const 
  {return apriori_covariance_;}

//-----------------------------------------------------------------------
/// Set apriori covariance value
//-----------------------------------------------------------------------

  void apriori_covariance(const blitz::Array<double, 2>& m) 
  {apriori_covariance_.reference(m.copy());}

//-----------------------------------------------------------------------
/// First guess value.
//-----------------------------------------------------------------------

  const blitz::Array<double, 1>& initial_guess() const 
  {
    if(initial_guess_.rows() > 0)
      return initial_guess_;
    else
      return apriori_;
  }

//-----------------------------------------------------------------------
/// Set first guess value.
//-----------------------------------------------------------------------

  void initial_guess(const blitz::Array<double, 1>& v) 
  {initial_guess_.reference(v.copy());}

  virtual void print(std::ostream& Os) const;

private:
  blitz::Array<double, 1> apriori_, initial_guess_;
  blitz::Array<double, 2> apriori_covariance_;
};
}
#endif
