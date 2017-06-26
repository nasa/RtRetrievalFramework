#ifndef STOKES_COEFFICIENT_CONSTANT_H
#define STOKES_COEFFICIENT_CONSTANT_H
#include "stokes_coefficient_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the stokes coefficient portion of the
  state. This particular implementation just uses constant values, 
  i.e. the state vector has no effect on the stokes coefficients.
*******************************************************************/
class StokesCoefficientConstant : public StokesCoefficientImpBase {
public:
  StokesCoefficientConstant(const blitz::Array<double, 2>& Stokes_coeff);
  virtual ~StokesCoefficientConstant() {}
  virtual void print(std::ostream& Os) const;
  virtual boost::shared_ptr<StokesCoefficient> clone() const;
  // Update the stokes coefficients.
  void set_stokes_coefficient(const blitz::Array<double, 2> Stokes_coeff)
  {
    if(Stokes_coeff.rows() != stokes_coeff.rows() ||
       Stokes_coeff.cols() != stokes_coeff.cols())
      throw Exception("Stokes_coeff is not the expected size");
    stokes_coeff = Stokes_coeff;
  }
protected:
  virtual void calc_stokes_coeff() const { 
    // Nothing to do
  }
};
}
#endif

