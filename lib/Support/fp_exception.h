#ifndef FP_EXCEPTION_H
#define FP_EXCEPTION_H
#include "printable.h"
#include <boost/backtrace.hpp>
#include <sstream>		// Definition of ostringstream.
#include <gsl/gsl_errno.h>

namespace FullPhysics {

/****************************************************************//**
  This is the base of the exception hierarchy for Full Physics code. 
  This can be written to like a stream to contain information about the
  exception. This is derived from the standard library std::exception
*******************************************************************/

class Exception: public std::exception, public Printable<Exception>,
		 public boost::backtrace {
public:
//-----------------------------------------------------------------------
/// Default constructor. Can give an optional string describing
/// the error.
//-----------------------------------------------------------------------

  Exception(const std::string& W = "") 
  { 
    // This reserve shouldn't really be necessary, but on a Mac
    // 10.4.11 using gcc 4.0.1, there is some kind of bug where we get
    // a "Double free" error when printing in Ruby. I never tracked
    // exactly where this occurred, but it was somewhere in the
    // iostream library when the buffer of os was resized. We just
    // reserve enough space up front so this isn't an issue. Since
    // this only gets called when an exception occurs, there shouldn't
    // be much of a performance issue with this.
    std::string buf("blah");
    buf.reserve(1000);
    s_.str(buf);
    s_ << W;  
  }

//-----------------------------------------------------------------------
/// Copy constructor.
//-----------------------------------------------------------------------

    Exception(const Exception& E)
    {
      try {
	s_ << E.s_.str();
      } catch(...) {		// Ignore all errors.
      }
    }

//-----------------------------------------------------------------------
/// Destructor.
//-----------------------------------------------------------------------
    
    virtual ~Exception() throw() {}

//-----------------------------------------------------------------------
/// Write to exception what() string.
//-----------------------------------------------------------------------
    
    template<class T> inline Exception& operator<<(const T& V)
    {
      s_ << V;
      return *this;
    }

//-----------------------------------------------------------------------
/// Print out description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  {
    Os << "Full Physics Exception:\n"
       << "=========================\n" 
       << what() << "\n"
       << "=========================\n" 
       << "Backtrace:\n"
       << "=========================\n" 
       << boost::trace(*this) << "\n"
       << "=========================\n";
  }

//-----------------------------------------------------------------------
/// Description of what the error is.
//-----------------------------------------------------------------------

    virtual const char* what() const throw()
    {

//-----------------------------------------------------------------------
// We can't just do s_.str().c_str(), because the temporary variable
// returned by str() disappears when we exit this function. Instead,
// we use a scratch variable that has the lifetime of this object.
//-----------------------------------------------------------------------
      try {
	scratch_ = s_.str();
      } catch(...) {		// If an error condition occurs,
				// ignore it and return a null
				// pointer.
	return 0;
      }
      return scratch_.c_str();
    }
  private:
    std::ostringstream s_;
    mutable std::string scratch_;
};


/** \defgroup Error Error checking routines */
/*@{*/
//-----------------------------------------------------------------------
/// Range check
//-----------------------------------------------------------------------

template <class T> inline void range_check_template(
const T&	   Val,		// Value to be checked.
const T&	   Min,		// Minimum allowed value.
const T&	   Max,		// Maximum allowed value.
const char*        File,
int                Line
)
{
  if(Val < Min ||
     !(Val < Max)) {
    Exception e;
    e << "Out of range error in file " << File << " at line " << Line << "\n"
      << "Value:           " << Val << "\n"
      << "Minimum allowed: " << Min << "\n"
      << "Maximum allowed: " << Max;
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Range check
//-----------------------------------------------------------------------

#define range_check(V, Min, Max) \
      FullPhysics::range_check_template(V, Min, Max, __FILE__, __LINE__)

//-----------------------------------------------------------------------
/// Range check
//-----------------------------------------------------------------------

template <class T> inline void range_min_check_template(
const T&	   Val,		// Value to be checked.
const T&	   Min,		// Minimum allowed value.
const char*        File,
int                Line
)
{
  if(Val < Min) {
    Exception e;
    e << "Out of range error in file " << File << " at line " << Line << "\n"
      << "Value:           " << Val << "\n"
      << "Minimum allowed: " << Min << "\n";
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Range check
//-----------------------------------------------------------------------

#define range_min_check(V, Min) \
      FullPhysics::range_min_check_template(V, Min, __FILE__, __LINE__)

//-----------------------------------------------------------------------
/// Range check
//-----------------------------------------------------------------------

template <class T> inline void range_max_check_template(
const T&	   Val,		// Value to be checked.
const T&	   Max,		// Maximum allowed value.
const char*        File,
int                Line
)
{
  if(!(Val < Max)) {
    Exception e;
    e << "Out of range error in file " << File << " at line " << Line << "\n"
      << "Value:           " << Val << "\n"
      << "Maximum allowed: " << Max << "\n";
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Range check
//-----------------------------------------------------------------------

#define range_max_check(V, Max) \
      FullPhysics::range_max_check_template(V, Max, __FILE__, __LINE__)

//-----------------------------------------------------------------------
/// Turn off gsl errors abort
//-----------------------------------------------------------------------

inline void no_gsl_abort()
{
  gsl_set_error_handler_off();
}

//-----------------------------------------------------------------------
/// Check for gsl errors
//-----------------------------------------------------------------------

inline void gsl_check_func
(int status,
const char* File,
int Line
)
{
  if(status != 0) {
    Exception e;
    e << "GSL error in file " << File << " at line " << Line << "\n"
      << "GSL error: " << gsl_strerror(status) << "\n";
    throw e;
  }
}

//-----------------------------------------------------------------------
/// GSL check
//-----------------------------------------------------------------------

#define gsl_check(status) \
      FullPhysics::gsl_check_func(status, __FILE__, __LINE__)

}
/*@}*/
#endif
