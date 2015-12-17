#include "fts_run_log.h"
#include "fp_exception.h"
#include <fstream>
#include <boost/regex.hpp>
using namespace boost::posix_time;
using namespace boost::gregorian;

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(FtsRunLog)
REGISTER_LUA_END()

REGISTER_LUA_CLASS(FtsRunLogRecord)
.def_readonly("outside_pressure", &FtsRunLogRecord::outside_pressure)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<FtsRunLogRecord>, FtsRunLogVector)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Read a FtsRunLogRecord from a HDF file group 
/// The datasets read match the same output by the L2 single sounding file
//-----------------------------------------------------------------------

FtsRunLog::FtsRunLog(const HdfFile& Hfile, const std::string& Group_name, const std::vector<std::string>& Band_names) {
  for(int sidx = 0; sidx < (int) Band_names.size(); sidx++) {
    FtsRunLogRecord rec;
    rec.spectrum_name = Band_names[sidx];

    TinyVector<int, 2> start(0, sidx);
    TinyVector<int, 2> size(1, 1);
    
    rec.time = Time::time_pgs( Hfile.read_field<double,2>(Group_name + "/time_pgs", start, size)(0,0) );
    rec.latitude = Hfile.read_field<double,2>(Group_name + "/latitude", start, size)(0,0);
    rec.longitude = Hfile.read_field<double,2>(Group_name + "/longitude", start, size)(0,0);
    rec.altitude = Hfile.read_field<double,2>(Group_name + "/altitude", start, size)(0,0);
    rec.solar_zenith = Hfile.read_field<double,2>(Group_name + "/solar_zenith", start, size)(0,0);
    rec.solar_azimuth = Hfile.read_field<double,2>(Group_name + "/solar_azimuth", start, size)(0,0);
    rec.zenith_offset = Hfile.read_field<double,2>(Group_name + "/zenith_offset", start, size)(0,0);
    rec.observer_sun_doppler_shift = Hfile.read_field<double,2>(Group_name + "/observer_sun_doppler_shift", start, size)(0,0);
    rec.optical_path_difference = Hfile.read_field<double,2>(Group_name + "/optical_path_difference", start, size)(0,0);
    rec.internal_fov = Hfile.read_field<double,2>(Group_name + "/internal_fov", start, size)(0,0);
    rec.external_fov = Hfile.read_field<double,2>(Group_name + "/external_fov", start, size)(0,0);
    rec.angular_misalignment = Hfile.read_field<double,2>(Group_name + "/angular_misalignment", start, size)(0,0);

    rec.zero_level_offset = Hfile.read_field<double,2>(Group_name + "/zero_level_offset", start, size)(0,0);
    rec.snr = Hfile.read_field<double,2>(Group_name + "/snr", start, size)(0,0);
    rec.inside_temperature = Hfile.read_field<double,2>(Group_name + "/inside_temperature", start, size)(0,0);
    rec.inside_pressure = Hfile.read_field<double,2>(Group_name + "/inside_pressure", start, size)(0,0);
    rec.inside_humidity = Hfile.read_field<double,2>(Group_name + "/inside_humidity", start, size)(0,0);
    rec.outside_temperature = Hfile.read_field<double,2>(Group_name + "/outside_temperature", start, size)(0,0);
    rec.outside_pressure = Hfile.read_field<double,2>(Group_name + "/outside_pressure", start, size)(0,0);
    rec.outside_humidity = Hfile.read_field<double,2>(Group_name + "/outside_humidity", start, size)(0,0);
    rec.solar_intensity_average = Hfile.read_field<double,2>(Group_name + "/solar_intensity_average", start, size)(0,0);
    rec.fractional_variation_solar_intensity = Hfile.read_field<double,2>(Group_name + "/fractional_variation_solar_intensity", start, size)(0,0);
    rec.wind_speed = Hfile.read_field<double,2>(Group_name + "/wind_speed", start, size)(0,0);
    rec.wind_direction = Hfile.read_field<double,2>(Group_name + "/wind_direction", start, size)(0,0);
    rec.laser_frequency = Hfile.read_field<double,2>(Group_name + "/laser_frequency", start, size)(0,0);
    rec.sun_tracker_frequency = Hfile.read_field<double,2>(Group_name + "/sun_tracker_frequency", start, size)(0,0);
    rec.airmass_independent_path_length = Hfile.read_field<double,2>(Group_name + "/airmass_independent_path_length", start, size)(0,0);
    rec.apodization_function = Hfile.read_field<std::string,2>(Group_name + "/apodization_function", start, size)(0,0);
    rec.spectrum_index = Hfile.read_field<std::string,2>(Group_name + "/spectrum_index", start, size)(0,0);

    run_log_record[rec.spectrum_name] = rec;
  }
}

//-----------------------------------------------------------------------
/// Read a FtsRunLogRecord from the given stream.
//-----------------------------------------------------------------------

std::istream& FullPhysics::operator>>(std::istream& is, FtsRunLogRecord& rec)
{
  is >> rec.spectrum_name;
  int year, doy;
  double tday;
  is >> year >> doy >> tday;
  int hr = (int) floor(tday);
  int mnt = (int) floor((tday - hr) * 60.0);
  int scnd = (int) floor(((tday - hr) * 60.0 - mnt) * 60.0);
  double frac_sec = ((tday - hr) * 60.0 - mnt) * 60.0 - scnd;
  ptime t(date(year, 1, 1) + date_duration(doy - 1), 
          hours(hr) + minutes(mnt) + seconds(scnd) +
          microseconds((long) floor(frac_sec * 1e6 + 0.5)));
  rec.time = Time(t);
  is >> rec.latitude  //oblat obs_info(num)%sounding_lat_topo,
     >> rec.longitude //oblon obs_info(num)%sounding_lon_topo,
     >> rec.altitude //obalt obs_info(num)%sounding_alt_topo,
     >> rec.solar_zenith //asza  obs_info(num)%solar_zenith_topo,
     >> rec.zenith_offset //zenoff zenoff
     >> rec.solar_azimuth //azim
     >> rec.observer_sun_doppler_shift //osds
     >> rec.optical_path_difference //opd state%instrument%ils(num)%coefs(1, 1)
     >> rec.internal_fov //fovi
     >> rec.external_fov //fovo
     >> rec.angular_misalignment //amal
     >> rec.index_first //ifirst fts_aux(num)%ifirst,
     >> rec.index_last //ilast fts_aux(num)%ilast
     >> rec.spacing_raw_spectrum //graw state%instrument%dispersion(2, num),
     >> rec.length_attached_header //possp
     >> rec.byte_per_word //bytepw
     >> rec.zero_level_offset //zoff
     >> rec.snr //snr
     >> rec.apodization_function //apf
     >> rec.inside_temperature //tins
     >> rec.inside_pressure //pins
     >> rec.inside_humidity //hins
     >> rec.outside_temperature
     >> rec.outside_pressure
     >> rec.outside_humidity
     >> rec.solar_intensity_average //sia
     >> rec.fractional_variation_solar_intensity //fvsi
     >> rec.wind_speed //wspd
     >> rec.wind_direction //wdir
     >> rec.laser_frequency //lasf
     >> rec.sun_tracker_frequency //wavtkr
     >> rec.airmass_independent_path_length //aipl
     >> rec.spectrum_index; //unqspcid

  // Try and obtain spectrum index from spectrum_name
  // This assumes that the convention of using the spectrum
  // index at the end of the filename
  if (rec.spectrum_index == "") {
    rec.spectrum_index = rec.spectrum_name.substr(rec.spectrum_name.size()-3,3);
  }

  return is;
}

//-----------------------------------------------------------------------
/// Read the given file.
//-----------------------------------------------------------------------

FtsRunLog::FtsRunLog(const std::string& Fname)
: fname(Fname)
{
  std::ifstream in(Fname.c_str());
  if(!in.good())
    throw Exception("Trouble reading FTS run log file " + fname);
  in.exceptions(std::ifstream::badbit);
  int line = 0;                        // Keep track of the line we are
                                       // reading so we can report it in any
                                       // error message
  std::string ln;

  // First line is header telling number of lines to skip in header and
  // number of columns
  int skip_header, number_columns;
  in >> skip_header >> number_columns;
  line++;
  
  // Skip lines to get to first data row
  for(int s = 0; s < skip_header; s++) {
    getline(in, ln);
    line++;
  }

  // Last line skipped should be names of columns
  // Count them so that it can be tested that the number of columns
  // match what we expect to parse
  boost::sregex_token_iterator icur(ln.begin(), ln.end(), boost::regex("\\s+"), -1);
  boost::sregex_token_iterator iend;
  int num_columns = 0;
  while(icur != iend) {
    num_columns++;
    icur++;
  }

  try {
    while(!in.eof() && in.good()) {
      getline(in, ln);
      line++;
      // If the first character in the line is ":", then skip the line.
      if(in.good() &&
         ln[0] != ':') {
        if(num_columns < 37)
          throw Exception("We only support runlog files with 37+ columns.");
        // Strip off a "-", "+", or " " in the first character.
        if(ln[0] == '-' || ln[0] == '+' || ln[0] == ' ')
          ln = ln.substr(1);
        // We don't currently handle the tab separated data.
        if(ln.find('\t') != std::string::npos)
          throw Exception("We don't currently support tab delimited data");
        std::istringstream ins(ln);
        FtsRunLogRecord rec;
        ins >> rec;
        run_log_record[rec.spectrum_name] = rec;
      }
    }
  } catch(Exception& exc) {
    exc << "\nWhile reading FTS run log File: " << fname << "\n"
        << "Line: " << line;
    throw;
  }
  
}

//-----------------------------------------------------------------------
/// Return the record for the given spectrum name.
//-----------------------------------------------------------------------

const FtsRunLogRecord& FtsRunLog::read(const std::string& spectrum_name) const
{
  std::map<std::string, FtsRunLogRecord>::const_iterator i = 
    run_log_record.find(spectrum_name);
  if(i == run_log_record.end()) {
    Exception e;
    e << "The spectrum '" << spectrum_name << "' was not found in the\n"
      << "the FTS run log file '" << fname << "'";
    throw e;
  }
  return i->second;
}
