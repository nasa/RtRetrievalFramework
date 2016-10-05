#ifndef ACCUMULATED_TIMER_H
#define ACCUMULATED_TIMER_H
#include "printable.h"
#include "logger.h"
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <boost/timer.hpp>

namespace FullPhysics {
class FunctionTimer;
class FunctionTimerR;

/****************************************************************//**
  This is a simple timer class that can be used to accumulate the time
  spent in multiple calls to a function. You get a function timer from
  this class, which keeps track of the time that object exists and
  adds it to the elapsed time.
*******************************************************************/
class AccumulatedTimer : public Printable<AccumulatedTimer> {
public:

//-----------------------------------------------------------------------
/// Create an AccumulatedTimer. The description is printed out when
/// you print this object.
//-----------------------------------------------------------------------

  AccumulatedTimer(const std::string& Desc) : elapsed_(0.0), desc(Desc) {}

//-----------------------------------------------------------------------
/// Total elapsed time.
//-----------------------------------------------------------------------
  
  double elapsed() const {return elapsed_; }
  
//-----------------------------------------------------------------------
/// Reset elapsed time to 0.
//-----------------------------------------------------------------------
  
  void reset_elapsed() {elapsed_ = 0.0;}

//-----------------------------------------------------------------------
/// Function timer
//-----------------------------------------------------------------------

  FunctionTimer function_timer(bool Auto_log = false) const;
  void print(std::ostream& Os) const 
  { Os << desc << " elapsed time " << elapsed(); }
private:
  mutable double elapsed_;
  std::string desc;
  friend class FunctionTimerR;
};

/****************************************************************//**
  Helper class for AccumulatedTimer
*******************************************************************/
class FunctionTimerR : boost::noncopyable {
public:
  FunctionTimerR(const AccumulatedTimer& At, bool Auto_log)  
    : at(At), auto_log(Auto_log) {}
  ~FunctionTimerR() 
  { 
    at.elapsed_ += t.elapsed(); 
    if(auto_log) {
      Logger::info() << "Current: " << at.desc << " elapsed time " 
		     << t.elapsed() << "\n";
      Logger::info() << "Total:   " << at << "\n";
    }
  }
private:
  const AccumulatedTimer& at;
  boost::timer t;
  bool auto_log;
};

/****************************************************************//**
  Helper class for AccumulatedTimer
*******************************************************************/

class FunctionTimer {
public:
  FunctionTimer(const AccumulatedTimer& At, bool Auto_log) : 
    p(new FunctionTimerR(At, Auto_log)) {}
private:
  boost::shared_ptr<FunctionTimerR> p;
};

inline FunctionTimer AccumulatedTimer::function_timer(bool Auto_log) const 
{ return FunctionTimer(*this, Auto_log); }
}
#endif
