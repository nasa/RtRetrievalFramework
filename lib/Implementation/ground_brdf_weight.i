// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ground_brdf_weight.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground)
%base_import(sub_state_vector_array)
%import "double_with_unit.i"
%import "array_with_unit.i"

%fp_shared_ptr(FullPhysics::GroundBrdfWeight);
namespace FullPhysics {
class GroundBrdfWeight: public SubStateVectorArray<Ground> {
public:
  GroundBrdfWeight(const blitz::Array<double, 2>& Spec_coeffs,
              const blitz::Array<bool,2>& Flag,
              const ArrayWithUnit<double, 1>& Ref_points,
              const std::vector<std::string>& Desc_band_names,
              const std::string& Weighted_brdf_name);

  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
  virtual const int number_spectrometer() const;
  virtual const int number_params() const;
  virtual const ArrayAd<double, 1> weight_coefficients(const int spec_index) const;
  virtual const blitz::Array<double, 2> weight_covariance(const int spec_index) const;
  virtual const DoubleWithUnit reference_point(const int spec_index) const;
  virtual boost::shared_ptr<Ground> clone() const;
  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;
  virtual std::string desc() const;
};
}
