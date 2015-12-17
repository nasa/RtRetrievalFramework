#ifndef OCO_ECMWF_H
#define OCO_ECMWF_H
#include "ecmwf.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements the OCO specific ECMWF reading 
  functionality.
*******************************************************************/

class OcoEcmwf : public Ecmwf {
public:
  OcoEcmwf(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& Hdf_sounding_id);
  ~OcoEcmwf() {}

  // Define how to read various items
  blitz::Array<double, 1> specific_humidity(const blitz::Array<double, 1>& Pressure_level) const
  { return read_and_interpolate("ECMWF/specific_humidity_profile_ecmwf", Pressure_level); }
 
  ArrayAd<double, 1> specific_humidity(const ArrayAd<double, 1>& Pressure_level) const
  { return read_and_interpolate("ECMWF/specific_humidity_profile_ecmwf", Pressure_level); }

  void temperature_grid(blitz::Array<double, 1>& Pressure, blitz::Array<double, 1>& T) const
  { read("ECMWF/temperature_profile_ecmwf", Pressure, T); }

  void specific_humidity_grid(blitz::Array<double, 1>& Pressure, blitz::Array<double, 1>& H) const
  { read("ECMWF/specific_humidity_profile_ecmwf", Pressure, H); }

  double surface_pressure() const
  { return read("ECMWF/surface_pressure_ecmwf"); }

  double windspeed() const
  { return sqrt( sqr(read("ECMWF/windspeed_u_ecmwf")) + sqr(read("ECMWF/windspeed_v_ecmwf")) ); }
  
  blitz::Array<double, 1> temperature(const blitz::Array<double, 1>& Pressure_level) const
  { return read_and_interpolate("ECMWF/temperature_profile_ecmwf", Pressure_level); }

  virtual ArrayAd<double, 1> temperature(const ArrayAd<double, 1>& Pressure_level) const
  { return read_and_interpolate("ECMWF/temperature_profile_ecmwf", Pressure_level); }

  void print(std::ostream& Os) const { Os << "OcoEcmwf"; }

private:

  //-----------------------------------------------------------------------
  /// OCO specific ECMWF reader routines
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
