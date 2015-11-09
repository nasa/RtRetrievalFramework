#ifndef ACOS_ECMWF_H
#define ACOS_ECMWF_H
#include "ecmwf.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements the ACOS specific ECMWF reading 
  functionality.
*******************************************************************/

class AcosEcmwf : public Ecmwf {
public:
  AcosEcmwf(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& 
	Hdf_sounding_id, bool Avg_sounding_number);
  AcosEcmwf(const std::string& Fname, const HeritageFile& Run_file);
  ~AcosEcmwf() {}

  // Define how to read various items
  blitz::Array<double, 1> specific_humidity(const blitz::Array<double, 1>& Pressure_level) const
  { return read_and_interpolate("ecmwf/specific_humidity", Pressure_level); }
 
  ArrayAd<double, 1> specific_humidity(const ArrayAd<double, 1>& Pressure_level) const
  { return read_and_interpolate("ecmwf/specific_humidity", Pressure_level); }

  void temperature_grid(blitz::Array<double, 1>& Pressure, blitz::Array<double, 1>& T) const
  { read("ecmwf/temperature", Pressure, T); }

  void specific_humidity_grid(blitz::Array<double, 1>& Pressure, blitz::Array<double, 1>& H) const
  { read("ecmwf/specific_humidity", Pressure, H); }

  double surface_pressure() const
  { return read("ecmwf/surface_pressure"); }

  double windspeed_u() const
  { return read("ecmwf/windspeed_u"); }

  double windspeed_v() const
  { return read("ecmwf/windspeed_v"); }
   
  blitz::Array<double, 1> temperature(const blitz::Array<double, 1>& Pressure_level) const
  { return read_and_interpolate("ecmwf/temperature", Pressure_level); }

  ArrayAd<double, 1> temperature(const ArrayAd<double, 1>& Pressure_level) const
  { return read_and_interpolate("ecmwf/temperature", Pressure_level); }

  void print(std::ostream& Os) const { Os << "AcosEcmwf"; }

private:

  //-----------------------------------------------------------------------
  /// ACOS specific ECMWF reader routines
  //-----------------------------------------------------------------------
  
  double read(const std::string& Field) const;
  void read(const std::string& Field, blitz::Array<double, 1>& P, 
	    blitz::Array<double, 1>& V) const;

  HdfFile h;
  boost::shared_ptr<HdfSoundingId> hsid;
  bool average_sounding_number;
};
}
#endif
