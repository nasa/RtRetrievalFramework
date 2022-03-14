#ifndef TEMPERATURE_LEVEL_SHAPE_H
#define TEMPERATURE_LEVEL_SHAPE_H

#include "temperature_imp_base.h"

namespace FullPhysics {

/****************************************************************//**
  This class maintains the temperature portion of the state. The
  temperature is retrieved as a set of scalings applied to a set of
  shape profiles combined with an initial guess temperature profile.

  The shape is computed as follows for each level i
  temp_final(i) = temp_base(i) + scaling1 * shape1(i) + scaling2 * shape2(i) + ... scalingN * shapeN(i)

  Where scalingN are the scale factors and shapeN are the shape
  profile values.
*******************************************************************/

class TemperatureLevelShape: public TemperatureImpBase {
public:
    TemperatureLevelShape(const blitz::Array<double, 1> Temp_base,
                          const blitz::Array<double, 2> Shape_profile,
                          const blitz::Array<double, 1> Shape_scaling,
                          const boost::shared_ptr<Pressure>& Press,
                          const blitz::Array<bool, 1>& Shape_flag);

    virtual ~TemperatureLevelShape() = default;

    virtual void print(std::ostream& Os) const;

    virtual boost::shared_ptr<Temperature> clone(const boost::shared_ptr<Pressure>& Press) const;

    virtual boost::shared_ptr<Temperature> clone() const { return clone(press->clone()); }

    virtual std::string sub_state_identifier() const
    {
        return "temperature_level_shape";
    }

    virtual std::string state_vector_name_i(int i) const;

    //-----------------------------------------------------------------------
    /// Base temperature profile associated with the pressure profile, values are in 
    /// Kelvin
    //-----------------------------------------------------------------------

    virtual blitz::Array<double, 1> temperature_base() const
    { return temp_base; }

    //-----------------------------------------------------------------------
    /// Shape profiles combined with base temperatue, dimenison are:
    /// N_level x N_scaling
    //-----------------------------------------------------------------------

    virtual blitz::Array<double, 2> shape_profile() const
    { return shape_prof; }

    //-----------------------------------------------------------------------
    /// Scale factors combined with the shape profiles with dimension:
    /// N_scaling 
    //-----------------------------------------------------------------------

    virtual blitz::Array<double, 1> shape_scaling() const
    { return coeff.value(); }

    //-----------------------------------------------------------------------
    /// Pressure levels that serve as the grid for the temperature values in
    /// units of Pascals
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> pressure_profile() const
    { return press->pressure_grid().value.value(); }

    virtual ArrayWithUnit<double, 1> important_pressure_level() const
    {
        return ArrayWithUnit<double, 1>(pressure_profile(), units::Pa);
    }
protected:
    void calc_temperature_grid() const;
private:
    blitz::Array<double, 1> temp_base;
    blitz::Array<double, 2> shape_prof;
};
}
#endif
