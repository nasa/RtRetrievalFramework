// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "dispersion_polynomial_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "dispersion_polynomial.i"
%fp_shared_ptr(FullPhysics::DispersionPolynomialOutput);

namespace FullPhysics {
class DispersionPolynomialOutput : public RegisterOutputBase {
public:
  DispersionPolynomialOutput(const boost::shared_ptr<DispersionPolynomial>& D,
			     const std::string& Hdf_band_name);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}


