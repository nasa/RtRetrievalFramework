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

class GroundBreon: public SubStateVectorArray<Ground> {
public:
    GroundBreon(const double Amplitude, const double Asymmetry, const double Geometric,
                const bool Ampl_flag, const bool Asym_flag, const bool Geom_flag);
    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
    virtual const AutoDerivative<double> overall_amplitude() const;
    virtual const AutoDerivative<double> asymmetry_parameter() const;
    virtual const AutoDerivative<double> geometric_factor() const;
    virtual const double refractive_index(const int Spec_idx) const;
    virtual const std::string breon_type() const = 0;
    virtual boost::shared_ptr<Ground> clone() const = 0;
    virtual std::string state_vector_name_i(int i) const;
    virtual void print(std::ostream& Os) const;
    virtual std::string desc() const;
};


class GroundBreonVeg: public GroundBreon {
public:
    GroundBreonVeg(const double Amplitude, const double Asymmetry, const double Geometric,
                const bool Ampl_flag, const bool Asym_flag, const bool Geom_flag);
    virtual const std::string breon_type() const;
    virtual boost::shared_ptr<Ground> clone() const;
};

class GroundBreonSoil: public GroundBreon {
public:
    GroundBreonSoil(const double Amplitude, const double Asymmetry, const double Geometric,
                const bool Ampl_flag, const bool Asym_flag, const bool Geom_flag);
    virtual const std::string breon_type() const;
    virtual boost::shared_ptr<Ground> clone() const;
};

}
