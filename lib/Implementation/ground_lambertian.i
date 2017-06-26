// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ground_lambertian.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground)
%base_import(sub_state_vector_array)
%import "double_with_unit.i"
%import "array_with_unit.i"

%fp_shared_ptr(FullPhysics::GroundLambertian);
namespace FullPhysics {
class GroundLambertian: public SubStateVectorArray<Ground> {
public:
  GroundLambertian(const blitz::Array<double, 2>& Spec_coeffs,
                   const blitz::Array<bool,2>& Flag, 
                   const ArrayWithUnit<double, 1>& Ref_points,
                   const std::vector<std::string>& Desc_band_names);

  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
  virtual const AutoDerivative<double> albedo(const DoubleWithUnit wave_point, const int spec_index) const;
  virtual const int number_spectrometer() const;
  virtual const int number_params() const;
  virtual const ArrayAd<double, 1> albedo_coefficients(const int spec_index) const;
  virtual const blitz::Array<double, 2> albedo_covariance(const int spec_index) const;
  virtual const DoubleWithUnit reference_point(const int spec_index) const;
  virtual boost::shared_ptr<Ground> clone() const;
  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;
  virtual std::string desc() const;
};
}
