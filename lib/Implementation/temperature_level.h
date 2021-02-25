#ifndef TEMPERATURE_LEVEL_H
#define TEMPERATURE_LEVEL_H

#include "temperature_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the temperature portion of the state. The
  temperature is retrieved through a set of levels related to pressures
  on the same grid.
*******************************************************************/
class TemperatureLevel: public TemperatureImpBase {
public:
    TemperatureLevel(const blitz::Array<double, 1> Temp,
                     const boost::shared_ptr<Pressure>& Press,
                     const blitz::Array<bool, 1>& Temp_flag);

    virtual ~TemperatureLevel() = default;

    virtual void print(std::ostream& Os) const;

    virtual boost::shared_ptr<Temperature> clone(const boost::shared_ptr<Pressure>& Press) const;

    virtual boost::shared_ptr<Temperature> clone() const { return clone(press->clone()); }

    virtual std::string sub_state_identifier() const
    {
        return "temperature_level";
    }

    virtual std::string state_vector_name_i(int i) const;

    //-----------------------------------------------------------------------
    /// Full temperature associate with the pressure profile, values are in 
    /// Kelvin
    //-----------------------------------------------------------------------

    virtual blitz::Array<double, 1> temperature_profile() const
    { return coeff.value(); }

    //-----------------------------------------------------------------------
    /// Pressure levels that serve as the grid for the temperature values in
    /// units of Pascals
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> pressure_profile() const
    { return pressure->pressure_grid().value.value(); }

    virtual ArrayWithUnit<double, 1> important_pressure_level() const
    {
        return ArrayWithUnit<double, 1>(pressure_profile(), units::Pa);
    }
protected:
    void calc_temperature_grid() const;
private:
    boost::shared_ptr<Pressure> pressure;
};
}
#endif
