#ifndef STOKES_COEFFICIENT_FRACTION_OUTPUT_H
#define STOKES_COEFFICIENT_FRACTION_OUTPUT_H
#include "register_output_base.h"
#include "stokes_coefficient_fraction.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the StokesCoefficientFraction class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the StokesCoefficientFraction class.
*******************************************************************/
class StokesCoefficientFractionOutput : public RegisterOutputBase {
public:
  StokesCoefficientFractionOutput(const boost::shared_ptr<StokesCoefficientFraction>& F, int Spec_index, const std::string& Hdf_band_name) 
    : f(F), spec_index(Spec_index), hdf_band_name(Hdf_band_name) {}
  virtual ~StokesCoefficientFractionOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<StokesCoefficientFraction> f;
  int spec_index;
  std::string hdf_band_name;
};
}
#endif
