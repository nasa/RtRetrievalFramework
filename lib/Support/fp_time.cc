#include "fp_time.h"
#include <cmath>

using namespace FullPhysics;
using namespace boost::posix_time;
using namespace boost::gregorian;

//-----------------------------------------------------------------------
/// Create from a boost posix time.
//-----------------------------------------------------------------------

Time::Time(const boost::posix_time::ptime& t)
{
  const ptime epoch(date(1970,1,1), time_duration(0,0,0,0));
  unix_time_ = (t - epoch).ticks() / 
    ((double) time_duration::ticks_per_second());
}

//-----------------------------------------------------------------------
/// Convert to a boost ptime.
//-----------------------------------------------------------------------

Time::operator boost::posix_time::ptime() const
{
  return ptime(date(1970,1,1), 
    time_duration(0,0,0,
     floor(unix_time() * time_duration::ticks_per_second() + 0.5)));
}


//-----------------------------------------------------------------------
/// Calculate the fractional day of the year. This is the number of
/// days, including the fractional part of the time of day, since
/// midnight January 1 of the year.
//-----------------------------------------------------------------------

double Time::frac_day_of_year() const
{
  ptime pt(*this);
  return pt.date().day_of_year() + pt.time_of_day().total_microseconds() /
    (24 * 60 * 60 * 1e6);
}

//-----------------------------------------------------------------------
/// Calculate the fractional year. This is the year plus the fractional
/// part of it.
//-----------------------------------------------------------------------

double Time::frac_year() const
{
  ptime pt(*this);
  // Fractional part of the year. Do it this way to account for leap years 
  // Use total_microseconds for greater precision
  double frac = double(( pt - ptime(date(pt.date().year(), 1, 1)) ).total_microseconds()) /
      double(( ptime(date(pt.date().year(), 12, 31)) - ptime(date(pt.date().year(), 1, 1)) ).total_microseconds());
  return pt.date().year() + frac;

}

//-----------------------------------------------------------------------
/// Parse CCSDS format time (e.g., "1996-07-03T04:13:57.987654Z")
//-----------------------------------------------------------------------

Time Time::parse_time(const std::string& Time_string)
{
  std::string tstring(Time_string);

  // boost::posix_time::time_from_string can't directly read CCSDS
  // format strings, but it can read something very close to
  // this. Massage the string into a format that time_from_string can
  // read.
  if(tstring.find('T') != std::string::npos)
    tstring[tstring.find('T')] = ' ';
  tstring = tstring.substr(0, tstring.find('Z'));
  return Time(boost::posix_time::time_from_string(tstring));
}

//-----------------------------------------------------------------------
/// Convert to CCSDS format.
//-----------------------------------------------------------------------

std::string Time::to_string() const
{
  std::string t = to_iso_extended_string(ptime(*this));
  // ptime doesn't write exactly CCSDS format. Massage this to make it
  // the right format.
  return t + "Z";
}
