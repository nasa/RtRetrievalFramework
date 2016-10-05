#ifndef DISPERSION_POLYNOMIAL_H
#define DISPERSION_POLYNOMIAL_H
#include "dispersion.h"
#include "unit.h"
#include "sub_state_vector_array.h"

namespace FullPhysics {
/****************************************************************//**
  This is an implementation of Dispersion that uses a polynomial
  expression to calculate the wavenumbers.

  Note that there are two minor variations of the dispersion
  polynomial. The first wavenumber returned can either be the
  polynomial evaluated at the value of "1", or a value of "0". By
  convention, the polynomial is 1 based for GOSAT and OCO, but 0 based
  for FTS.
*******************************************************************/
class DispersionPolynomial: public SubStateVectorArray<Dispersion> {
public:
  DispersionPolynomial(const blitz::Array<double, 1>& Coeff, 
		       const blitz::Array<bool, 1>& Used_flag,
		       const Unit& Coeff_unit,
		       const std::string& Band_name,
		       int Number_pixel,
		       bool Is_one_based);
  DispersionPolynomial(const blitz::Array<double, 1>& Coeff, 
		       const blitz::Array<bool, 1>& Used_flag,
		       const std::string& Coeff_unit_name,
		       const std::string& Band_name,
		       int Number_pixel, 
		       bool Is_one_based);
  virtual ~DispersionPolynomial() {}

//-----------------------------------------------------------------------
/// Dispersion offset. This is just coeff(0), but we wrap this for use
/// by DispersionPolynomialOutput
//-----------------------------------------------------------------------

  double dispersion_offset() const { return coeff.value()(0); }

//-----------------------------------------------------------------------
/// Dispersion spacing. This is just coeff(1), but we wrap this for use
/// by DispersionPolynomialOutput
//-----------------------------------------------------------------------

  double dispersion_spacing() const { return coeff.value()(1); }

//-----------------------------------------------------------------------
/// Dispersion offset uncertainty. This is just sqrt(Cov(0,0)), but we
/// wrap this for use  by DispersionPolynomialOutput
//-----------------------------------------------------------------------

  double dispersion_offset_uncertainty() const 
  { 
    if(sv_cov_sub.rows() < 1)
      return 0;
    double t = sv_cov_sub(0,0);
    return (t < 0 ? 0 : sqrt(t)); 
  }

//-----------------------------------------------------------------------
/// Dispersion spacing uncertainty. This is just sqrt(Cov(1,1)), but we
/// wrap this for use  by DispersionPolynomialOutput
//-----------------------------------------------------------------------

  double dispersion_spacing_uncertainty() const 
  { 
    if(sv_cov_sub.rows() < 2)
      return 0;
    double t = sv_cov_sub(1,1);
    return (t < 0 ? 0 : sqrt(t)); 
  }

  virtual boost::shared_ptr<Dispersion> clone() const;
  virtual std::string state_vector_name_i(int i) const;
  virtual SpectralDomain pixel_grid() const;
  virtual void print(std::ostream& Os) const;
private:
  void initialize();
  bool is_one_based;
  Unit coeff_unit;
  std::string band_name_;
  /// This is an array like 1,2,3 ... number_pixel. This is used by 
  /// pixel_grid
  blitz::Array<double, 1> index_array;
  // Very similar to index_array, but always 1 based.
  blitz::Array<int, 1> spectral_index;
};
}
#endif
