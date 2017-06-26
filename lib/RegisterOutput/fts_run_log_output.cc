#include "fts_run_log_output.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(FtsRunLogOutput, RegisterOutputBase)
.def(luabind::constructor<const std::vector<FtsRunLogRecord>&>())
REGISTER_LUA_END()
#endif

// See base class for description.
void FtsRunLogOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  Array<std::string, 1> spectrum_index((int) run_log.size());
  Array<double, 1> time_pgs(spectrum_index.shape());
  Array<std::string, 1> time_string(spectrum_index.shape());
  Array<double, 1> latitude(spectrum_index.shape());
  Array<double, 1> longitude(spectrum_index.shape());
  Array<double, 1> altitude(spectrum_index.shape());
  Array<double, 1> solar_zenith(spectrum_index.shape());
  Array<double, 1> solar_azimuth(spectrum_index.shape());
  Array<double, 1> zenith_offset(spectrum_index.shape());
  Array<double, 1> observer_sun_doppler_shift(spectrum_index.shape());
  Array<double, 1> optical_path_difference(spectrum_index.shape());
  Array<double, 1> internal_fov(spectrum_index.shape());
  Array<double, 1> external_fov(spectrum_index.shape());
  Array<double, 1> angular_misalignment(spectrum_index.shape());
  Array<double, 1> zero_level_offset(spectrum_index.shape());
  Array<double, 1> snr(spectrum_index.shape());  
  Array<double, 1> inside_temperature(spectrum_index.shape());
  Array<double, 1> inside_pressure(spectrum_index.shape());
  Array<double, 1> inside_humidity(spectrum_index.shape());
  Array<double, 1> outside_temperature(spectrum_index.shape());
  Array<double, 1> outside_pressure(spectrum_index.shape());
  Array<double, 1> outside_humidity(spectrum_index.shape());
  Array<double, 1> solar_intensity_average(spectrum_index.shape());
  Array<double, 1> fractional_variation_solar_intensity(spectrum_index.shape());
  Array<double, 1> wind_speed(spectrum_index.shape());
  Array<double, 1> wind_direction(spectrum_index.shape());
  Array<double, 1> laser_frequency(spectrum_index.shape());
  Array<double, 1> sun_tracker_frequency(spectrum_index.shape());
  Array<double, 1> airmass_independent_path_length(spectrum_index.shape());
  Array<std::string, 1> apodization_function(spectrum_index.shape());

  for(int i = 0; i < spectrum_index.rows(); ++i) {
    spectrum_index(i) = run_log[i].spectrum_index;
    time_pgs(i) = run_log[i].time.pgs_time();
    time_string(i) = run_log[i].time.to_string();
    latitude(i) = run_log[i].latitude;
    longitude(i) = run_log[i].longitude;
    altitude(i) = run_log[i].altitude;
    solar_zenith(i) = run_log[i].solar_zenith;
    solar_azimuth(i) = run_log[i].solar_azimuth;
    zenith_offset(i) = run_log[i].zenith_offset;
    observer_sun_doppler_shift(i) = run_log[i].observer_sun_doppler_shift;
    optical_path_difference(i) = run_log[i].optical_path_difference;
    internal_fov(i) = run_log[i].internal_fov;
    external_fov(i) = run_log[i].external_fov;
    angular_misalignment(i) = run_log[i].angular_misalignment;
    zero_level_offset(i) = run_log[i].zero_level_offset;
    snr(i) = run_log[i].snr;
    inside_temperature(i) = run_log[i].inside_temperature;
    inside_pressure(i) = run_log[i].inside_pressure;
    inside_humidity(i) = run_log[i].inside_humidity;
    outside_temperature(i) = run_log[i].outside_temperature;
    outside_pressure(i) = run_log[i].outside_pressure;
    outside_humidity(i) = run_log[i].outside_humidity;
    solar_intensity_average(i) = run_log[i].solar_intensity_average;
    fractional_variation_solar_intensity(i) = 
      run_log[i].fractional_variation_solar_intensity;
    wind_speed(i) = run_log[i].wind_speed;
    wind_direction(i) = run_log[i].wind_direction;
    laser_frequency(i) = run_log[i].laser_frequency;
    sun_tracker_frequency(i) = run_log[i].sun_tracker_frequency;
    airmass_independent_path_length(i) = 
      run_log[i].airmass_independent_path_length;
    apodization_function(i) = run_log[i].apodization_function;
  }
  out->register_data_source("FtsRunLog/spectrum_index", spectrum_index);
  out->register_data_source("FtsRunLog/time_pgs", time_pgs);
  out->register_data_source("FtsRunLog/time_string", time_string);
  out->register_data_source("FtsRunLog/latitude", latitude);
  out->register_data_source("FtsRunLog/longitude", longitude);
  out->register_data_source("FtsRunLog/altitude", altitude);
  out->register_data_source("FtsRunLog/solar_zenith", solar_zenith);
  out->register_data_source("FtsRunLog/solar_azimuth", solar_azimuth);
  out->register_data_source("FtsRunLog/zenith_offset", zenith_offset);
  out->register_data_source("FtsRunLog/observer_sun_doppler_shift",   
                           observer_sun_doppler_shift);
  out->register_data_source("FtsRunLog/optical_path_difference",
                           optical_path_difference);
  out->register_data_source("FtsRunLog/internal_fov", internal_fov);
  out->register_data_source("FtsRunLog/external_fov", external_fov);
  out->register_data_source("FtsRunLog/angular_misalignment",
                           angular_misalignment);
  out->register_data_source("FtsRunLog/zero_level_offset", zero_level_offset);
  out->register_data_source("FtsRunLog/snr", snr);
  out->register_data_source("FtsRunLog/inside_temperature", inside_temperature);
  out->register_data_source("FtsRunLog/inside_pressure", inside_pressure);
  out->register_data_source("FtsRunLog/inside_humidity", inside_humidity);
  out->register_data_source("FtsRunLog/outside_temperature",
                           outside_temperature);
  out->register_data_source("FtsRunLog/outside_pressure", outside_pressure);
  out->register_data_source("FtsRunLog/outside_humidity", outside_humidity);
  out->register_data_source("FtsRunLog/solar_intensity_average",
                           solar_intensity_average);
  out->register_data_source("FtsRunLog/fractional_variation_solar_intensity",
                           fractional_variation_solar_intensity);
  out->register_data_source("FtsRunLog/wind_speed", wind_speed);
  out->register_data_source("FtsRunLog/wind_direction", wind_direction);
  out->register_data_source("FtsRunLog/laser_frequency", laser_frequency);
  out->register_data_source("FtsRunLog/sun_tracker_frequency", 
                           sun_tracker_frequency);
  out->register_data_source("FtsRunLog/airmass_independent_path_length",
                           airmass_independent_path_length);
  out->register_data_source("FtsRunLog/apodization_function",
                           apodization_function);
}
