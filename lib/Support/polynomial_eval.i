// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "polynomial_eval.h"
%}
%base_import(generic_object)
%import "array_ad.i"
%import "auto_derivative.i"

%fp_shared_ptr(FullPhysics::Poly1d)

namespace FullPhysics {
class Poly1d : public GenericObject {
public:
  Poly1d(const ArrayAd<double, 1>& Coefficients);
  AutoDerivative<double> operator()(const AutoDerivative<double>& Value) const;
  ArrayAd<double, 1> operator()(const ArrayAd<double, 1>& Array) const;
  std::string print_to_string() const;
};
}
