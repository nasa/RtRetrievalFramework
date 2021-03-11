// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "ils_gaussian.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}
%base_import(ils_function)
%import "array_ad.i"
%import "auto_derivative.i"
%fp_shared_ptr(FullPhysics::IlsGaussian);

namespace FullPhysics {
class IlsGaussian : public IlsFunction {
public:
  IlsGaussian(double a, const std::string& Band_name, 
	      const std::string& Hdf_band_name);
  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn, ArrayAd<double, 1>& OUTPUT,
   bool jac_optimization=false) const;
};
}
