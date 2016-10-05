
#ifndef FTS_RUN_LOG_H
#define FTS_RUN_LOG_H
#include "fp_time.h"
#include "hdf_file.h"

namespace FullPhysics {
/****************************************************************//**
  This is a single FTS run log record.
*******************************************************************/
struct FtsRunLogRecord : public Printable<FtsRunLogRecord> {
  std::string spectrum_name;
  Time time;
  /// Observation latitude in degrees.
  double latitude;
  /// Observation longitude in degrees.
  double longitude;
  /// Observation altitude, in degrees.
  double altitude;
  /// Astronomical solar zenith angle (unrefracted)
  double solar_zenith;
  /// Zenith angle pointing offset (degree)
  double zenith_offset;
  /// Solar azimuth angle
  double solar_azimuth;
  /// Observer-Sun Doppler Shift (ppm)
  double observer_sun_doppler_shift;
  /// Optical path difference (cm) of interferogram
  double optical_path_difference;
  /// Internal angular diameter of FOV (radians)
  double internal_fov;
  /// External angular diameter of FOV (radians)
  double external_fov;
  /// Angular misalignment of interferometer (radians)
  double angular_misalignment;
  /// Index of first spectral point in disk file
  int index_first;
  /// Index of last spectral point in disk file
  int index_last;
  // Spacing of raw spectrum (cm-1) from GETINFO
  double spacing_raw_spectrum;
  /// Length of attached header in bytes
  int length_attached_header;
  /// Bytes per data word (usually 2 or 4).
  int byte_per_word;
  /// Zero level offset (dimensionless fraction)
  double zero_level_offset;
  /// Signal-to-Noise Ratio (dimensionless)
  double snr;
  /// Apodization function;
  std::string apodization_function;
  /// Inside temperature
  double inside_temperature;
  /// Inside pressure
  double inside_pressure;
  /// Inside humidity
  double inside_humidity;
  /// Outside temperature
  double outside_temperature;
  /// Outside pressure;
  double outside_pressure;
  /// Outside humidity;
  double outside_humidity;
  /// Solar Intensity (Average)
  double solar_intensity_average;
  /// Fractional Variation in Solar Intensity
  double fractional_variation_solar_intensity;
  /// Wind Speed (m/s)
  double wind_speed;
  /// Wind Direction (deg)
  double wind_direction;
  /// Laser Frequency (e.g. 15798 cm-1)
  double laser_frequency;
  /// Suntracker frequency (active tracking)
  double sun_tracker_frequency;
  /// Airmass-Independent Path Length (km)
  double airmass_independent_path_length;
  /// Spectrum index.
  std::string spectrum_index;

  /// Initializes all the records to a sane default
  FtsRunLogRecord() {
    spectrum_name = "";
    time = Time();
    latitude = 0.0;
    longitude = 0.0;
    altitude = 0.0;
    solar_zenith = 0.0;
    zenith_offset = 0.0;
    solar_azimuth = 0.0;
    observer_sun_doppler_shift = 0.0;
    optical_path_difference = 0.0;
    internal_fov = 0.0;
    external_fov = 0.0;
    angular_misalignment = 0.0;
    index_first = 0;
    index_last = 0;
    spacing_raw_spectrum = 0.0;
    length_attached_header = 0;
    byte_per_word = 0;
    zero_level_offset = 0.0;
    snr = 0.0;
    apodization_function = "";
    inside_temperature = 0.0;
    inside_pressure = 0.0;
    inside_humidity = 0.0;
    outside_temperature = 0.0;
    outside_pressure = 0.0;
    outside_humidity = 0.0;
    solar_intensity_average = 0.0;
    fractional_variation_solar_intensity = 0.0;
    wind_speed = 0.0;
    wind_direction = 0.0;
    laser_frequency = 0.0;
    sun_tracker_frequency = 0.0;
    airmass_independent_path_length = 0.0;
    spectrum_index = "";
  }

  virtual void print(std::ostream& Os) const {Os << "FtsRunLogRecord";}  
};

std::istream& operator>>(std::istream& is, FtsRunLogRecord& rec);

/****************************************************************//**
  This reads a FTS run log file. This is just a simple text file, this
  class handles reading this data.

  There are two different types of run logs, either space
  delimited or tab delimited. We currently only support the space
  delimited (support for tab delimited would really just require
  adding unit test data, it isn't all that different from the space
  delimited).

  There are also several formats for the space delimited files, which
  vary depending on the size of the a line. We only read the latest
  "New GDS-format", which has a line length of 300.
*******************************************************************/
  class FtsRunLog : public Printable<FtsRunLog> {
public:
  FtsRunLog(const std::string& Fname);
  FtsRunLog(const HdfFile& Hfile, const std::string& Group_name, const std::vector<std::string>& Band_names);
  const FtsRunLogRecord& read(const std::string& spectrum_name) const;
  virtual void print(std::ostream& Os) const {Os << "FtsRunLog";}
private:
  std::map<std::string, FtsRunLogRecord> run_log_record;
  std::string fname;
};

}
#endif
