#ifndef OCO_SIM_MET_ECMWF_H
#define OCO_SIM_MET_ECMWF_H
#include "ecmwf.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements the OCO specific ECMWF reading 
  functionality.

  This reads the simulator meteorology files. This are similar to 
  the OCO ECMWF, but the fields have different names.

  Note that the actual OCO simulator used "scene" files, which are
  somewhat lime the meteorology but in a different format, and with
  a different number of levels (not the normal 91 ECMWF). The
  meteorology files are this scene information resampled to the 91
  levels. 
*******************************************************************/

class OcoSimMetEcmwf : public Ecmwf {
public:
  OcoSimMetEcmwf(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& Hdf_sounding_id);
  ~OcoSimMetEcmwf() {}

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

  double windspeed() const
  { return sqrt( sqr(read("ecmwf/windspeed_u")) + sqr(read("ecmwf/windspeed_v")) ); }
  
  blitz::Array<double, 1> temperature(const blitz::Array<double, 1>& Pressure_level) const
  { return read_and_interpolate("ecmwf/temperature", Pressure_level); }

  virtual ArrayAd<double, 1> temperature(const ArrayAd<double, 1>& Pressure_level) const
  { return read_and_interpolate("ecmwf/temperature", Pressure_level); }

  void print(std::ostream& Os) const { Os << "OcoSimMetEcmwf"; }

private:

  //-----------------------------------------------------------------------
  /// OCO simulator meteorology specific ECMWF reader routines
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
