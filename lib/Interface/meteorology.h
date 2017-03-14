#ifndef METEOROLOGY_H
#define METEOROLOGY_H
#include "hdf_sounding_id.h"
#include "hdf_file.h"
#include "array_ad.h"
#include "printable.h"

namespace FullPhysics {
/****************************************************************//**
 Defines the interface for supplying meteorological data. Routines
 are provided to interpolate the values to a different pressure
 levels grid for values with a profile.
*******************************************************************/

class Meteorology : public Printable<Meteorology> {
public:
    virtual ~Meteorology() {}

    //-----------------------------------------------------------------------
    /// Pressure levels in Pascals used for provided meteorological data
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> pressure_levels() const = 0;

    //-----------------------------------------------------------------------
    /// Specific humidty on the meteorological pressure levels
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> specific_humidity() const = 0;

    //-----------------------------------------------------------------------
    /// Specific humidity interpolated to the requested pressure levels
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> specific_humidity(const blitz::Array<double, 1>& Pressure_level) const;

    //-----------------------------------------------------------------------
    /// Volume mixing ratio for a particular species on the meteorological 
    /// pressure levels. Species name should use the standard molecular compound
    /// naming convention. Ie, oxygen = O2, water = H2O, etc..
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> vmr(const std::string& Species) const = 0;

    //-----------------------------------------------------------------------
    /// Volume mixing ratio for a particular species interpolated to the 
    /// requested pressure levels
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> vmr(const std::string& Species, const blitz::Array<double, 1>& Pressure_level) const;

    //-----------------------------------------------------------------------
    /// Temperature profile in Kelvins on the meteorological pressure levels
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> temperature() const = 0;

    //-----------------------------------------------------------------------
    /// Temperature profile in Kelvins interpolated to the requested pressure levels
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> temperature(const blitz::Array<double, 1>& Pressure_level) const;

    //-----------------------------------------------------------------------
    /// Surface pressure in Pascals
    //-----------------------------------------------------------------------
    virtual double surface_pressure() const = 0;

    //-----------------------------------------------------------------------
    /// Windspeed in m/s for the surface
    //-----------------------------------------------------------------------
    virtual double windspeed() const = 0;

    void print(std::ostream& Os) const { Os << "Meteorology"; }

protected:

    //-----------------------------------------------------------------------
    /// Interpolates a profile of data from the internal pressure grid
    /// to the supplied one. Uses log-log interpolation.
    //-----------------------------------------------------------------------
    blitz::Array<double, 1> interpolate_to_grid(const blitz::Array<double, 1>& Profile, const blitz::Array<double, 1>& Dest_pressure_levels) const;

};
}
#endif
