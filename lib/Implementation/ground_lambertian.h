#ifndef GROUND_LAMBERTIAN_H
#define GROUND_LAMBERTIAN_H

#include "ground.h"
#include "array_with_unit.h"
#include "sub_state_vector_array.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a Lambertian albedo as a ground type. 
*******************************************************************/
class GroundLambertian: public SubStateVectorArray<Ground> {

public:

  GroundLambertian(const blitz::Array<double, 2>& Spec_coeffs,
                   const blitz::Array<bool,2>& Flag, 
                   const ArrayWithUnit<double, 1>& Ref_points,
                   const std::vector<std::string>& Desc_band_names);

  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

  virtual const AutoDerivative<double> albedo(const DoubleWithUnit wave_point, const int spec_index) const;

  virtual const int number_spectrometer() const { return desc_band_names.size(); }
  virtual const int number_params() const { return coefficient().value().rows() / number_spectrometer(); }

  virtual const ArrayAd<double, 1> albedo_coefficients(const int spec_index) const;
  virtual const blitz::Array<double, 2> albedo_covariance(const int spec_index) const;

  /// Center wavelength that spectrally dependent parameter is referenced to
  virtual const DoubleWithUnit reference_point(const int spec_index) const { return reference_points(spec_index); }

  virtual boost::shared_ptr<Ground> clone() const;

  virtual std::string state_vector_name_i(int i) const;

  virtual void print(std::ostream& Os) const;

  virtual std::string desc() const { return "GroundLambertian"; }

protected:

  GroundLambertian(const blitz::Array<double, 1>& Spec_coeffs,
                   const blitz::Array<bool, 1>& Flag, 
                   const ArrayWithUnit<double, 1>& Ref_points,
                   const std::vector<std::string>& Desc_band_names);

private:
  
  ArrayWithUnit<double, 1> reference_points;
  std::vector<std::string> desc_band_names;

};
}
#endif
