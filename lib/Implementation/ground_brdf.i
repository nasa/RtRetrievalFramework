// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ground_brdf.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground)
%base_import(sub_state_vector_array)

%fp_shared_ptr(FullPhysics::GroundBrdfVeg);
%fp_shared_ptr(FullPhysics::GroundBrdfSoil);

namespace FullPhysics {

class GroundBrdfVeg: public Ground {
public:
    GroundBrdfVeg(const blitz::Array<double, 2>& Coeffs,
                const blitz::Array<bool, 2>& Flag,
                const std::vector<std::string>& Desc_band_names);
    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
    virtual const int number_spectrometer() const { return desc_band_names.size(); }
    virtual const AutoDerivative<double> rahman_factor(const int spec_index) const;
    virtual const AutoDerivative<double> overall_amplitude(const int spec_index) const;
    virtual const AutoDerivative<double> asymmetry_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> geometric_factor(const int spec_index) const;
    virtual const AutoDerivative<double> breon_factor(const int spec_index) const;
    virtual void rahman_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void overall_amplitude(const int spec_index, const AutoDerivative<double>& val);
    virtual void asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void geometric_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void breon_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual const double refractive_index(const int Spec_idx) const;
    virtual const std::string breon_type() const;
    virtual boost::shared_ptr<Ground> clone() const;
    virtual std::string state_vector_name_i(int i) const;
    virtual void print(std::ostream& Os) const;
    virtual std::string desc() const;
};

class GroundBrdfSoil: public Ground {
public:
    GroundBrdfSoil(const blitz::Array<double, 2>& Coeffs,
                const blitz::Array<bool, 2>& Flag,
                const std::vector<std::string>& Desc_band_names);
    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
    virtual const int number_spectrometer() const { return desc_band_names.size(); }
    virtual const AutoDerivative<double> rahman_factor(const int spec_index) const;
    virtual const AutoDerivative<double> overall_amplitude(const int spec_index) const;
    virtual const AutoDerivative<double> asymmetry_parameter(const int spec_index) const;
    virtual const AutoDerivative<double> geometric_factor(const int spec_index) const;
    virtual const AutoDerivative<double> breon_factor(const int spec_index) const;
    virtual void rahman_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void overall_amplitude(const int spec_index, const AutoDerivative<double>& val);
    virtual void asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val);
    virtual void geometric_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual void breon_factor(const int spec_index, const AutoDerivative<double>& val);
    virtual const double refractive_index(const int Spec_idx) const;
    virtual const std::string breon_type() const;
    virtual boost::shared_ptr<Ground> clone() const;
    virtual std::string state_vector_name_i(int i) const;
    virtual void print(std::ostream& Os) const;
    virtual std::string desc() const;
};

}
