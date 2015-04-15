// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ground_breon.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground)
%base_import(sub_state_vector_array)

%fp_shared_ptr(FullPhysics::GroundBreonVeg);
%fp_shared_ptr(FullPhysics::GroundBreonSoil);

namespace FullPhysics {

class GroundBreonVeg: public Ground {
public:
    GroundBreonVeg(const double Amplitude, const double Asymmetry, const double Geometric,
                const bool Ampl_flag, const bool Asym_flag, const bool Geom_flag);
    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
    %python_attribute_with_set(overall_amplitude, AutoDerivative<double>);
    %python_attribute_with_set(asymmetry_parameter, AutoDerivative<double>);
    %python_attribute_with_set(geometric_factor, AutoDerivative<double>);
    virtual const double refractive_index(const int Spec_idx) const;
    virtual const std::string breon_type() const;
    virtual boost::shared_ptr<Ground> clone() const;
    virtual std::string state_vector_name_i(int i) const;
    virtual void print(std::ostream& Os) const;
    virtual std::string desc() const;
};

class GroundBreonSoil: public Ground {
public:
    GroundBreonSoil(const double Amplitude, const double Asymmetry, const double Geometric,
                const bool Ampl_flag, const bool Asym_flag, const bool Geom_flag);
    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
    %python_attribute_with_set(overall_amplitude, AutoDerivative<double>);
    %python_attribute_with_set(asymmetry_parameter, AutoDerivative<double>);
    %python_attribute_with_set(geometric_factor, AutoDerivative<double>);
    virtual const double refractive_index(const int Spec_idx) const;
    virtual const std::string breon_type() const;
    virtual boost::shared_ptr<Ground> clone() const;
    virtual std::string state_vector_name_i(int i) const;
    virtual void print(std::ostream& Os) const;
    virtual std::string desc() const;
};

}
