// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "fts_run_log.h"
%}
%base_import(generic_object)
%import "fp_time.i"
%import "hdf_file.i"
%fp_shared_ptr(FullPhysics::FtsRunLogRecord)
%fp_shared_ptr(FullPhysics::FtsRunLog)

namespace FullPhysics {
class FtsRunLogRecord : public GenericObject {
public:
  std::string spectrum_name;
  Time time;
  double latitude;
  double longitude;
  double altitude;
  double solar_zenith;
  double zenith_offset;
  double solar_azimuth;
  double observer_sun_doppler_shift;
  double optical_path_difference;
  double internal_fov;
  double external_fov;
  double angular_misalignment;
  int index_first;
  int index_last;
  double spacing_raw_spectrum;
  int length_attached_header;
  int byte_per_word;
  double zero_level_offset;
  double snr;
  std::string apodization_function;
  double inside_temperature;
  double inside_pressure;
  double inside_humidity;
  double outside_temperature;
  double outside_pressure;
  double outside_humidity;
  double solar_intensity_average;
  double fractional_variation_solar_intensity;
  double wind_speed;
  double wind_direction;
  double laser_frequency;
  double sun_tracker_frequency;
  double airmass_independent_path_length;
  std::string spectrum_index;
};

class FtsRunLog : public GenericObject {
public:
  FtsRunLog(const HdfFile& Hfile, const std::string& Group_name, const std::vector<std::string>& Band_names);
  FtsRunLog(const std::string& Fname);
  const FtsRunLogRecord& read(const std::string& spectrum_name) const;
};
}

