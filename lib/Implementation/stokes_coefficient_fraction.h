#ifndef STOKES_COEFFICIENT_FRACTION_H
#define STOKES_COEFFICIENT_FRACTION_H
#include "stokes_coefficient_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the stokes coefficient portion of the
  state. This particular implementation uses constant values, combined
  with fitting for the parallel polarization fraction.
*******************************************************************/
class StokesCoefficientFraction : public StokesCoefficientImpBase {
public:
  StokesCoefficientFraction(const blitz::Array<double, 2>& Stokes_coeff_parallel,
			    const blitz::Array<double, 1>& Coeffs,
			    const blitz::Array<bool, 1>& Flag);
  virtual ~StokesCoefficientFraction() {}
  virtual void print(std::ostream& Os) const;
  virtual boost::shared_ptr<StokesCoefficient> clone() const;
  virtual std::string state_vector_name_i(int i) const
  { return "Parallel Polarization Band " + 
      boost::lexical_cast<std::string>(i + 1); }

//-----------------------------------------------------------------------
/// f uncertainty. This is just sqrt(Cov(i,i)), but we
/// wrap this for use  by StokesCoefficientFractionOutput
//-----------------------------------------------------------------------

  double f_uncertainty(int Spec_index) const 
  { 
    if(sv_cov_sub.rows() < 1)
      return 0;
    double t = sv_cov_sub(Spec_index,Spec_index);
    return (t < 0 ? 0 : sqrt(t)); 
  }
  
protected:
  virtual void calc_stokes_coeff() const;
private:
  blitz::Array<double, 2> stokes_coeff_parallel;
};
}
#endif

