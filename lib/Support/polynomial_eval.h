#ifndef POLYNOMIAL_EVAL_H
#define POLYNOMIAL_EVAL_H

#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
 A one-dimensional polynomial class.

 A convenience class, used to encapsulate "natural" operations on
 polynomials so that said operations may take on their customary 
 form in code.

 Evaluation is done Horner's Scheme to reduce problems due to
 round off error and overflows.

 Inspired by numpy.poly1d

 Additionally handles ArrayAd data correctly.
*******************************************************************/
class Poly1d : public Printable<Poly1d> {
public:

  //-----------------------------------------------------------------------
  /// The polynomial's coefficients, in decreasing powers.
  //-----------------------------------------------------------------------
  Poly1d(const ArrayAd<double, 1>& Coefficients, const bool Decreasing_order = true) : coeffs_(Coefficients), decreasing_order_(Decreasing_order) {}

  //-----------------------------------------------------------------------
  /// Evaluate polynomial for a value
  //-----------------------------------------------------------------------
  double operator()(double Value) const;
  AutoDerivative<double> operator()(const AutoDerivative<double>& Value) const;

  //-----------------------------------------------------------------------
  /// Evaluate polynomial for an array
  //-----------------------------------------------------------------------
  blitz::Array<double,1> operator()(const blitz::Array<double,1>& Arr) const;
  ArrayAd<double, 1> operator()(const ArrayAd<double, 1>& Arr) const;
  
  virtual void print(std::ostream& Os) const;

private:

  blitz::Array<int,1> eval_indexes() const;

  ArrayAd<double, 1> coeffs_;
  bool decreasing_order_;
};

}
#endif
