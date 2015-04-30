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
    GroundBreon(const blitz::Array<double, 2>& Coeffs,
                const blitz::Array<bool, 2>& Flag,
                const std::vector<std::string>& Desc_band_names);

    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

    virtual const int number_spectrometer() const { return desc_band_names.size(); }

    // Rahman parameters
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
   
    /// Returns hard coded value of 1.5 since that is the value hardcoded into LIDORT
    virtual const double refractive_index(const int Spec_idx) const { return 1.5; }
  
    /// String describing which type of Breon surface type, also makes this class abstract
    virtual const std::string breon_type() const = 0;
    
    virtual boost::shared_ptr<Ground> clone() const = 0;
  
    virtual std::string state_vector_name_i(int i) const;
  
    virtual void print(std::ostream& Os) const;
  
    virtual std::string desc() const { return "GroundBreon"; }

protected:

    GroundBreon(const blitz::Array<double, 1>& Spec_coeffs,
                const blitz::Array<bool, 1>& Flag, 
                const std::vector<std::string>& Desc_band_names);

    std::vector<std::string> desc_band_names;
};


class GroundBreonVeg: public GroundBreon {
public:
    GroundBreonVeg(const blitz::Array<double, 2>& Coeffs,
                const blitz::Array<bool, 2>& Flag,
                const std::vector<std::string>& Desc_band_names) :
        GroundBreon(Coeffs, Flag, Desc_band_names) {}

    virtual const std::string breon_type() const { return "Vegetative"; }

    virtual boost::shared_ptr<Ground> clone() const {
      return boost::shared_ptr<Ground>(new GroundBreonVeg(coefficient().value(), used_flag_value(), desc_band_names));
    }
private:
    GroundBreonVeg(const blitz::Array<double, 1>& Spec_coeffs,
                const blitz::Array<bool, 1>& Flag, 
                const std::vector<std::string>& Desc_band_names) :
        GroundBreon(Spec_coeffs, Flag, Desc_band_names) {}
};

class GroundBreonSoil: public GroundBreon {
public:
    GroundBreonSoil(const blitz::Array<double, 2>& Coeffs,
                const blitz::Array<bool, 2>& Flag,
                const std::vector<std::string>& Desc_band_names) :
        GroundBreon(Coeffs, Flag, Desc_band_names) {}

    virtual const std::string breon_type() const { return "Soil"; }

    virtual boost::shared_ptr<Ground> clone() const {
      return boost::shared_ptr<Ground>(new GroundBreonSoil(coefficient().value(), used_flag_value(), desc_band_names));
    }
private:
    GroundBreonSoil(const blitz::Array<double, 1>& Spec_coeffs,
                const blitz::Array<bool, 1>& Flag, 
                const std::vector<std::string>& Desc_band_names) :
        GroundBreon(Spec_coeffs, Flag, Desc_band_names) {}
};

} // End of namespace

 #endif
