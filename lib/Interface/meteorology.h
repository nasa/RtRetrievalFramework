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
    /// 
    /// Functions for computing the VMR values are registered into
    /// vmr_func_map. This base class includes an implementation for H2O
    /// VMR from specific_humidity.
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> vmr(const std::string& Species) const;

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
    /// Windspeed magnitude in m/s for the surface
    //-----------------------------------------------------------------------
    virtual double windspeed() const;

    //-----------------------------------------------------------------------
    /// The U component windspeed in m/s
    //-----------------------------------------------------------------------
    virtual double windspeed_u() const = 0;

    //-----------------------------------------------------------------------
    /// The V component windspeed in m/s
    //-----------------------------------------------------------------------
    virtual double windspeed_v() const = 0;

    void print(std::ostream& Os) const { Os << "Meteorology"; }

protected:

    //-----------------------------------------------------------------------
    /// Return the H20 VMR. This is the specific_humidity converted to a
    /// volume mixing ratio.
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> h2o_vmr() const;

    //-----------------------------------------------------------------------
    /// Interpolates a profile of data from the internal pressure grid
    /// to the supplied one. Uses log-log interpolation.
    //-----------------------------------------------------------------------
    blitz::Array<double, 1> interpolate_to_grid(const blitz::Array<double, 1>& Profile, const blitz::Array<double, 1>& Dest_pressure_levels) const;

};
}
#endif
