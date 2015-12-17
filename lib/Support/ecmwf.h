#ifndef ECMWF_H
#define ECMWF_H
#include "hdf_sounding_id.h"
#include "hdf_file.h"
#include "array_ad.h"
#include "printable.h"

namespace FullPhysics {
/****************************************************************//**
  This class is used to read some of the fields from the ECMWF file,
  which can then be used for things such as the apriori.

  Since resampled ECMWF files can differ between instrument types,
  the read routines are pure virtual and need to be implemented 
  for the specifics of the instrument specific ECMWF files.
*******************************************************************/

class Ecmwf : public Printable<Ecmwf> {
public:
  virtual ~Ecmwf() {}

//-----------------------------------------------------------------------
/// Get the specific humidity at the given pressure levels. The ECMWF data
/// is on its own set of pressure levels, we interpolate to get this
/// at the requested levels.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> specific_humidity(const blitz::Array<double, 1>& Pressure_level) const = 0;
  virtual ArrayAd<double, 1> specific_humidity(const ArrayAd<double, 1>& Pressure_level) const = 0;

  blitz::Array<double, 1> h2o_vmr(const blitz::Array<double, 1>& Pressure_level) const;
  ArrayAd<double, 1> h2o_vmr(const ArrayAd<double, 1>& Pressure_level) const;

//-----------------------------------------------------------------------
/// Temperature grid on the ECMWF pressure grid. The temperature is in
/// Kelvin, and the Pressure is in pascals.
//-----------------------------------------------------------------------

  virtual void temperature_grid(blitz::Array<double, 1>& Pressure, blitz::Array<double, 1>& T) const = 0;

//-----------------------------------------------------------------------
/// Humidity on the ECMWF pressure grid. 
/// Pressure is in pascals.
//-----------------------------------------------------------------------

  virtual void specific_humidity_grid(blitz::Array<double, 1>& Pressure, blitz::Array<double, 1>& H) const = 0;

//-----------------------------------------------------------------------
/// Get the surface pressure from the Ecmwf file.
//-----------------------------------------------------------------------

  virtual double surface_pressure() const = 0;

//-----------------------------------------------------------------------
/// Get the windspeed from the Ecmwf file.
//-----------------------------------------------------------------------

  virtual double windspeed() const = 0;

//-----------------------------------------------------------------------
/// Get the temperature at the given pressure levels. The ECMWF data
/// is on its own set of pressure levels, we interpolate to get this
/// at the requested levels.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> temperature(const blitz::Array<double, 1>& Pressure_level) const = 0;
  virtual ArrayAd<double, 1> temperature(const ArrayAd<double, 1>& Pressure_level) const = 0;

  void print(std::ostream& Os) const { Os << "Ecmwf"; }

protected:
  //-----------------------------------------------------------------------
  /// Reader routines that need to be implemented for instrument specific
  /// ECMWF files
  //-----------------------------------------------------------------------
  
  virtual double read(const std::string& Field) const = 0;
  virtual void read(const std::string& Field, blitz::Array<double, 1>& P, 
		    blitz::Array<double, 1>& V) const = 0;

  //-----------------------------------------------------------------------
  /// Defines pressure based interpolation routines from ECMWF file data
  //-----------------------------------------------------------------------

  blitz::Array<double, 1> read_and_interpolate(const std::string& Field,
       const blitz::Array<double, 1>& Pressure_level) const;
  ArrayAd<double, 1> read_and_interpolate(const std::string& Field,
       const ArrayAd<double, 1>& Pressure_level) const;
};
}
#endif
