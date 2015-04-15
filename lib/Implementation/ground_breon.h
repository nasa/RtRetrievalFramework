#ifndef GROUND_BREON_H
#define GROUND_BREON_H

#include "ground.h"
#include "auto_derivative.h"
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
class GroundBreon: public SubStateVectorArray<Ground> {
public:
    GroundBreon(const double Amplitude, const double Asymmetry, const double Geometric,
                const bool Ampl_flag, const bool Asym_flag, const bool Geom_flag);

    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

    // Rahman parameters
    virtual const AutoDerivative<double> overall_amplitude() const;
    virtual const AutoDerivative<double> asymmetry_parameter() const;
    virtual const AutoDerivative<double> geometric_factor() const;

    virtual void overall_amplitude(const AutoDerivative<double>& val);
    virtual void asymmetry_parameter(const AutoDerivative<double>& val);
    virtual void geometric_factor(const AutoDerivative<double>& val);
   
    /// Returns hard coded value of 1.5 since that is the value hardcoded into LIDORT
    virtual const double refractive_index(const int Spec_idx) const { return 1.5; }
  
    /// String describing which type of Breon surface type, also makes this class abstract
    virtual const std::string breon_type() const = 0;
    
    virtual boost::shared_ptr<Ground> clone() const = 0;
  
    virtual std::string state_vector_name_i(int i) const;
  
    virtual void print(std::ostream& Os) const;
  
    virtual std::string desc() const { return "GroundBreon"; }

protected:
    std::vector<std::string> desc_band_names;
};


class GroundBreonVeg: public GroundBreon {
public:
    GroundBreonVeg(const double Amplitude, const double Asymmetry, const double Geometric,
                const bool Ampl_flag, const bool Asym_flag, const bool Geom_flag) :
        GroundBreon(Amplitude, Asymmetry, Geometric, Ampl_flag, Asym_flag, Geom_flag) {}

    virtual const std::string breon_type() const { return "Vegetative"; }

    virtual boost::shared_ptr<Ground> clone() const {
      return boost::shared_ptr<GroundBreonVeg>(new GroundBreonVeg(coefficient().value()(0), coefficient().value()(1), coefficient().value()(2),
                  used_flag_value()(0), used_flag_value()(1), used_flag_value()(2)));
    }
};

class GroundBreonSoil: public GroundBreon {
public:
    GroundBreonSoil(const double Amplitude, const double Asymmetry, const double Geometric,
                const bool Ampl_flag, const bool Asym_flag, const bool Geom_flag) :
        GroundBreon(Amplitude, Asymmetry, Geometric, Ampl_flag, Asym_flag, Geom_flag) {}

    virtual const std::string breon_type() const { return "Soil"; }

    virtual boost::shared_ptr<Ground> clone() const {
      return boost::shared_ptr<GroundBreonVeg>(new GroundBreonVeg(coefficient().value()(0), coefficient().value()(1), coefficient().value()(2),
                  used_flag_value()(0), used_flag_value()(1), used_flag_value()(2)));
    }
};

} // End of namespace

 #endif
