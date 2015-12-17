// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "stokes_coefficient_fraction.h"
%}

%base_import(stokes_coefficient_imp_base)
%fp_shared_ptr(FullPhysics::StokesCoefficientFraction);

namespace FullPhysics {
class StokesCoefficientFraction : public StokesCoefficientImpBase {
public:
  StokesCoefficientFraction(const blitz::Array<double, 2>& Stokes_coeff_parallel,
			    const blitz::Array<double, 1>& Coeffs,
			    const blitz::Array<bool, 1>& Flag);
  virtual boost::shared_ptr<StokesCoefficient> clone() const;
protected:
  virtual void calc_stokes_coeff() const;
};
}
