// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "stokes_coefficient_fraction_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "stokes_coefficient_fraction.i"
%fp_shared_ptr(FullPhysics::StokesCoefficientFractionOutput);

namespace FullPhysics {
class StokesCoefficientFractionOutput : public RegisterOutputBase {
public:
  StokesCoefficientFractionOutput
    (const boost::shared_ptr<StokesCoefficientFraction>& F,
     int Spec_index,
     const std::string& Hdf_band_name);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}


