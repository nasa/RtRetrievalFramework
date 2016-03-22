#ifndef FP_TIME_H
#define FP_TIME_H
#include "printable.h"
#include <boost/operators.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace FullPhysics {
/****************************************************************//**
  This is a simple time class.

  There are a few reasonable choices for expressing time information. 
  We could use TAI, GPS, the PGS toolkit. Each of these time system
  can be related to the other by a constant, since the only difference
  is the Epoch that time is measure against.

  For simplicity, we just use the standard unix time in seconds from
  January 1, 1970.

  This class abstracts out the representation we use for time. We
  supply conversions to the specific time systems for use in the cases
  that a specific system is needed (e.g., calling a PGS toolkit
  routine).

  We also supply methods for converting to and from a string
  representation of the time. We support the standard CCSDS format 
  (e.g., "1996-07-03T04:13:57.987654Z").
*******************************************************************/

class Time : public Printable<Time>,
             private boost::less_than_comparable<Time>,
	     private boost::addable<Time, double>,
	     private boost::subtractable<Time, double> {
public:

//-----------------------------------------------------------------------
/// Default constructor.
//-----------------------------------------------------------------------

  Time() { unix_time_ = 0; }

  Time(const boost::posix_time::ptime& t);
  operator boost::posix_time::ptime() const;

//-----------------------------------------------------------------------
/// Return time from given Unix time (epoch of 1970-01-01).
//-----------------------------------------------------------------------
  
  static Time time_unix(double unix_time) 
  {Time res; res.unix_time_ = unix_time; return res;}

//-----------------------------------------------------------------------
/// Return time from given PGS toolkit time (epoch of 1993-01-01).
//-----------------------------------------------------------------------
  
  static Time time_pgs(double pgs) 
  {Time res; res.unix_time_ = pgs + 725846400.0; return res;}

//-----------------------------------------------------------------------
/// Add given number of seconds to Time.
//-----------------------------------------------------------------------

  Time& operator+=(double T) {unix_time_ += T; return *this;}

//-----------------------------------------------------------------------
/// Subtract given number of seconds to Time.
//-----------------------------------------------------------------------

  Time& operator-=(double T) {unix_time_ -= T; return *this;}

//-----------------------------------------------------------------------
/// Give time in unix time, as a double (epoch 1970-01-01)
//-----------------------------------------------------------------------

  double unix_time() const { return unix_time_;}

//-----------------------------------------------------------------------
/// Give time in PGS toolkit time, as a double (epoch 1993-01-01)
//-----------------------------------------------------------------------

  double pgs_time() const { return unix_time_ - 725846400.0;}

  double frac_day_of_year() const;
  double frac_year() const;

  static Time parse_time(const std::string& Time_string);
  std::string to_string() const;

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

  void print(std::ostream& os) const { os << to_string(); }

private:
  double unix_time_;
};

//-----------------------------------------------------------------------
/// \ingroup Miscellaneous
///
/// Subtract two Times, giving the interval between them in seconds.
//-----------------------------------------------------------------------

inline double operator-(const Time& T1, const Time& T2) 
{ return T1.unix_time() - T2.unix_time(); }

//-----------------------------------------------------------------------
/// \ingroup Miscellaneous
/// Compare Times
///
/// We define <=, >= and > in terms of this operator.
//-----------------------------------------------------------------------

inline bool operator<(const Time& T1, const Time& T2)
{ return T1.unix_time() <  T2.unix_time(); }

}
#endif
