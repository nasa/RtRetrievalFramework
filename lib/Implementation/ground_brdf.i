// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ground_brdf.h"
#include "array_with_unit.h"
#include "double_with_unit.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground)
%base_import(array_with_unit)
%base_import(double_with_unit)
%base_import(sub_state_vector_array)

%fp_shared_ptr(FullPhysics::GroundBrdfVeg);
%fp_shared_ptr(FullPhysics::GroundBrdfSoil);

namespace FullPhysics {

class GroundBrdfVeg: public Ground {
public:
    GroundBrdfVeg(const blitz::Array<double, 2>& Coeffs,
                const blitz::Array<bool, 2>& Flag,
                const ArrayWithUnit<double, 1>& Ref_points, 
                const std::vector<std::string>& Desc_band_names);
    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
    virtual const int number_spectrometer() const { return desc_band_names.size(); }
    virtual const AutoDerivative<double> weight(const double wn, const int spec_index) const;
    virtual const AutoDerivative<double> weight_intercept(const int spec_index) const;
    virtual const AutoDerivative<double> weight_slope(const int spec_index) const;
  AutoDerivative<double> GroundBrdf::weight_coeff(const int spec_index, const int weight_index) const;
  ArrayAd<double, 1> weight_parameters(const int spec_index) const;
  virtual const AutoDerivative<double> rahman_factor(const int spec_index) const;
    virtual const AutoDerivative<double> hotspot_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> asymmetry_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> anisotropy_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> breon_factor(const int spec_index) const;
    virtual void weight_intercept(const int spec_index, const AutoDerivative<double>& val);
    virtual void weight_slope(const int spec_index, const AutoDerivative<double>& val);
    virtual void rahman_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void hotspot_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void anisotropy_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void breon_factor(const int spec_index, const AutoDerivative<double>& val);
    const blitz::Array<double, 2> brdf_covariance(const int spec_index) const;
    virtual const double refractive_index(const int Spec_idx) const;
    virtual const double black_sky_albedo(const int Spec_index, const double Sza);
    virtual const double kernel_value(const int Spec_index, const double Sza, const double Vza, const double Azm);
    virtual const std::string breon_type() const;
    virtual const DoubleWithUnit reference_point(const int spec_index) const;
    virtual boost::shared_ptr<Ground> clone() const;
    virtual std::string state_vector_name_i(int i) const;
    virtual void print(std::ostream& Os) const;
    virtual std::string desc() const;
};

class GroundBrdfSoil: public Ground {
public:
    GroundBrdfSoil(const blitz::Array<double, 2>& Coeffs,
                const blitz::Array<bool, 2>& Flag,
                const ArrayWithUnit<double, 1>& Ref_points,
                const std::vector<std::string>& Desc_band_names);
    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
    virtual const int number_spectrometer() const { return desc_band_names.size(); }
    virtual const AutoDerivative<double> weight(const double wn, const int spec_index) const;
    virtual const AutoDerivative<double> weight_intercept(const int spec_index) const;
    virtual const AutoDerivative<double> weight_slope(const int spec_index) const;
    virtual const AutoDerivative<double> rahman_factor(const int spec_index) const;
    virtual const AutoDerivative<double> hotspot_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> asymmetry_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> anisotropy_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> breon_factor(const int spec_index) const;
    virtual void weight_intercept(const int spec_index, const AutoDerivative<double>& val);
    virtual void weight_slope(const int spec_index, const AutoDerivative<double>& val);
    virtual void rahman_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void hotspot_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void anisotropy_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void breon_factor(const int spec_index, const AutoDerivative<double>& val);
    const blitz::Array<double, 2> brdf_covariance(const int spec_index) const;
    virtual const double refractive_index(const int Spec_idx) const;
    virtual const double black_sky_albedo(const int Spec_index, const double Sza);
    virtual const double kernel_value(const int Spec_index, const double Sza, const double Vza, const double Azm);
    virtual const std::string breon_type() const;
    virtual const DoubleWithUnit reference_point(const int spec_index) const;
    virtual boost::shared_ptr<Ground> clone() const;
    virtual std::string state_vector_name_i(int i) const;
    virtual void print(std::ostream& Os) const;
    virtual std::string desc() const;
};

}
