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
    virtual const AutoDerivative<double> rahman_factor(const int spec_index) const;
    virtual const AutoDerivative<double> overall_amplitude(const double wn, const int spec_index) const;
    virtual const AutoDerivative<double> overall_amplitude_intercept(const int spec_index) const;
    virtual const AutoDerivative<double> overall_amplitude_slope(const int spec_index) const;
    virtual const AutoDerivative<double> asymmetry_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> geometric_factor(const int spec_index) const;
    virtual const AutoDerivative<double> breon_factor(const int spec_index) const;
    virtual void rahman_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void overall_amplitude_intercept(const int spec_index, const AutoDerivative<double>& val);
    virtual void overall_amplitude_slope(const int spec_index, const AutoDerivative<double>& val);
    virtual void asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void geometric_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void breon_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual const double refractive_index(const int Spec_idx) const;
    virtual const double black_sky_albedo(const int Spec_index, const double Sza);
    virtual const double albedo(const int Spec_index, const double Sza, const double Vza, const double Azm, const blitz::Array<double, 1>& Stokes_coef);
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
    virtual const AutoDerivative<double> rahman_factor(const int spec_index) const;
    virtual const AutoDerivative<double> overall_amplitude(const double wn, const int spec_index) const;
    virtual const AutoDerivative<double> overall_amplitude_intercept(const int spec_index) const;
    virtual const AutoDerivative<double> overall_amplitude_slope(const int spec_index) const;
    virtual const AutoDerivative<double> asymmetry_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> geometric_factor(const int spec_index) const;
    virtual const AutoDerivative<double> breon_factor(const int spec_index) const;
    virtual void rahman_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void overall_amplitude_intercept(const int spec_index, const AutoDerivative<double>& val);
    virtual void overall_amplitude_slope(const int spec_index, const AutoDerivative<double>& val);
    virtual void asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void geometric_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void breon_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual const double refractive_index(const int Spec_idx) const;
    virtual const double black_sky_albedo(const int Spec_index, const double Sza);
    virtual const double albedo(const int Spec_index, const double Sza, const double Vza, const double Azm, const blitz::Array<double, 1>& Stokes_coef);
    virtual const std::string breon_type() const;
    virtual const DoubleWithUnit reference_point(const int spec_index) const;
    virtual boost::shared_ptr<Ground> clone() const;
    virtual std::string state_vector_name_i(int i) const;
    virtual void print(std::ostream& Os) const;
    virtual std::string desc() const;
};

}
