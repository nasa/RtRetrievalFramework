// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ils_function.h"
%}
%base_import(generic_object)

%import "auto_derivative.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::IlsFunction);

namespace FullPhysics {
class IlsFunction : public GenericObject {
public:
  virtual ~IlsFunction();
  std::string print_to_string() const;
  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn, ArrayAd<double, 1>& OUTPUT) const = 0;
  %python_attribute(band_name, virtual std::string);
  %python_attribute(hdf_band_name, virtual std::string);
};
}
