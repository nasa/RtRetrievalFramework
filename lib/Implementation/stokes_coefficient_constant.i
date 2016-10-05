// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "stokes_coefficient_constant.h"
%}

%base_import(stokes_coefficient_imp_base)
%fp_shared_ptr(FullPhysics::StokesCoefficientConstant);

namespace FullPhysics {
class StokesCoefficientConstant : public StokesCoefficientImpBase {
public:
  StokesCoefficientConstant(const blitz::Array<double, 2>& Stokes_coeff);
  virtual boost::shared_ptr<StokesCoefficient> clone() const;
  void set_stokes_coefficient(const blitz::Array<double, 2> Stokes_coeff);
protected:
  virtual void calc_stokes_coeff() const;
};
}
