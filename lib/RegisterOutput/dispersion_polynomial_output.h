#ifndef DISPERSION_POLYNOMIAL_OUTPUT_H
#define DISPERSION_POLYNOMIAL_OUTPUT_H
#include "register_output_base.h"
#include "dispersion_polynomial.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the DispersionPolynomial class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the DispersionPolynomial class.
*******************************************************************/
class DispersionPolynomialOutput : public RegisterOutputBase {
public:
  DispersionPolynomialOutput(const boost::shared_ptr<DispersionPolynomial>& D,
			     const std::string& Hdf_band_name) 
    : d(D), 
      hdf_band_name(Hdf_band_name) {}
  virtual ~DispersionPolynomialOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<DispersionPolynomial> d;
  std::string hdf_band_name;
};
}
#endif
