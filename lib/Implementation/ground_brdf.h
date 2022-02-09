#ifndef GROUND_BREON_H
#define GROUND_BREON_H

#include "ground.h"
#include "auto_derivative.h"
#include "array_with_unit.h"
#include "sub_state_vector_array.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements the Breon + Rahman ground type.
  The retrieved parameters are for the Rahman type. Refractive
  index is hardcoded at 1.5 as it is hard coded inside of LIDORT
  and l_rad.

  Base class for both Vegetative and Soil types which share the
  same mechanisms but are implemented differently in RT codes
*******************************************************************/
class GroundBrdf: public SubStateVectorArray<Ground> {
public:
    enum ParamIndex {
        RAHMAN_KERNEL_FACTOR_INDEX = 0,
        RAHMAN_OVERALL_AMPLITUDE_INDEX = 1,
        RAHMAN_ASYMMETRY_FACTOR_INDEX = 2,
        RAHMAN_GEOMETRIC_FACTOR_INDEX = 3,
        BREON_KERNEL_FACTOR_INDEX = 4,
        BRDF_WEIGHT_INTERCEPT_INDEX = 5,
        BRDF_WEIGHT_SLOPE_INDEX = 6
    };

    GroundBrdf(const blitz::Array<double, 2>& Coeffs,
               const blitz::Array<bool, 2>& Flag,
               const ArrayWithUnit<double, 1>& Ref_points,
               const std::vector<std::string>& Desc_band_names);

    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

    virtual const int number_spectrometer() const { return desc_band_names.size(); }

    virtual const int number_weight_parameters() const { return num_weight_params; }

    virtual const AutoDerivative<double> weight(const double wn, const int spec_index) const;

    virtual const AutoDerivative<double> rahman_factor(const int spec_index) const;
    virtual const AutoDerivative<double> hotspot_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> asymmetry_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> anisotropy_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> breon_factor(const int spec_index) const;
    virtual const AutoDerivative<double> weight_intercept(const int spec_index) const;
    virtual const AutoDerivative<double> weight_slope(const int spec_index) const;
    AutoDerivative<double> weight_coeff(const int spec_index, const int weight_index) const;
    ArrayAd<double, 1> weight_parameters(const int spec_index) const;

    virtual void rahman_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void hotspot_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void anisotropy_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void breon_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void weight_intercept(const int spec_index, const AutoDerivative<double>& val);
    virtual void weight_slope(const int spec_index, const AutoDerivative<double>& val);
    void weight_coeff(const int spec_index, const AutoDerivative<double>& val, const int weight_index);
    void weight_parameters(const int spec_index, const ArrayAd<double, 1>& val);

    const blitz::Array<double, 2> brdf_covariance(const int spec_index) const;

    /// Returns hard coded value of 1.5 since that is the value hardcoded into LIDORT
    virtual const double refractive_index(const int Spec_idx) const { return 1.5; }

    // Uses LIDORT to compute the black sky albedo from the parameters
    virtual const double black_sky_albedo(const int Spec_index, const double Sza) = 0;

    // Computes kernel value using parameters and specified geometry
    virtual const double kernel_value(const int Spec_index, const double Sza, const double Vza, const double Azm) = 0;

    /// String describing which type of Breon surface type, also makes this class abstract
    virtual const std::string breon_type() const = 0;

    /// Center wavelength that spectrally dependent parameter is referenced to
    virtual const DoubleWithUnit reference_point(const int spec_index) const { return reference_points(spec_index); }

    virtual boost::shared_ptr<Ground> clone() const = 0;

    virtual std::string state_vector_name_i(int i) const;

    virtual void print(std::ostream& Os) const;

    virtual std::string desc() const { return "GroundBrdf"; }

protected:
    GroundBrdf(const blitz::Array<double, 1>& Spec_coeffs,
               const blitz::Array<bool, 1>& Flag,
               const ArrayWithUnit<double, 1>& Ref_points,
               const std::vector<std::string>& Desc_band_names);

    int num_weight_params;
    int num_coeff;

    ArrayWithUnit<double, 1> reference_points;
    std::vector<std::string> desc_band_names;

    // Helper function for routines that call fortran codes
    blitz::Array<double, 1> black_sky_params(const int Spec_index);
    blitz::Array<double, 1> kernel_value_params(const int Spec_index);
};


class GroundBrdfVeg: public GroundBrdf {
public:
    GroundBrdfVeg(const blitz::Array<double, 2>& Coeffs,
                  const blitz::Array<bool, 2>& Flag,
                  const ArrayWithUnit<double, 1>& Ref_points,
                  const std::vector<std::string>& Desc_band_names) :
        GroundBrdf(Coeffs, Flag, Ref_points, Desc_band_names) {}

    virtual const double black_sky_albedo(const int Spec_index, const double Sza);
    virtual const double kernel_value(const int Spec_index, const double Sza, const double Vza, const double Azm);
    virtual const std::string breon_type() const { return "Vegetative"; }

    virtual boost::shared_ptr<Ground> clone() const {
      return boost::shared_ptr<Ground>(new GroundBrdfVeg(coefficient().value(), used_flag_value(), reference_points, desc_band_names));
    }
private:
    GroundBrdfVeg(const blitz::Array<double, 1>& Spec_coeffs,
                  const blitz::Array<bool, 1>& Flag,
                  const ArrayWithUnit<double, 1>& Ref_points,
                  const std::vector<std::string>& Desc_band_names) :
        GroundBrdf(Spec_coeffs, Flag, Ref_points, Desc_band_names) {}
};

class GroundBrdfSoil: public GroundBrdf {
public:
    GroundBrdfSoil(const blitz::Array<double, 2>& Coeffs,
                   const blitz::Array<bool, 2>& Flag,
                   const ArrayWithUnit<double, 1>& Ref_points,
                   const std::vector<std::string>& Desc_band_names) :
        GroundBrdf(Coeffs, Flag, Ref_points, Desc_band_names) {}

    virtual const double black_sky_albedo(const int Spec_index, const double Sza);
    virtual const double kernel_value(const int Spec_index, const double Sza, const double Vza, const double Azm);
    virtual const std::string breon_type() const { return "Soil"; }

    virtual boost::shared_ptr<Ground> clone() const {
      return boost::shared_ptr<Ground>(new GroundBrdfSoil(coefficient().value(), used_flag_value(), reference_points, desc_band_names));
    }
private:
    GroundBrdfSoil(const blitz::Array<double, 1>& Spec_coeffs,
                   const blitz::Array<bool, 1>& Flag,
                   const ArrayWithUnit<double, 1>& Ref_points,
                   const std::vector<std::string>& Desc_band_names) :
        GroundBrdf(Spec_coeffs, Flag, Ref_points, Desc_band_names) {}
};

} // End of namespace

 #endif
