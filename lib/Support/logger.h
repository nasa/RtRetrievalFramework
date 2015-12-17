#ifndef LOGGER_H
#define LOGGER_H
#include "printable.h"
#include <boost/shared_ptr.hpp>
#include <sstream>

namespace FullPhysics {
/****************************************************************//**
  The actual implementation of the Logger.
*******************************************************************/
class LogImp : public Printable<LogImp> {
public:
  virtual ~LogImp() {}
  enum log_level {DEBUG = 4, INFO=3, WARNING=2, ERROR=1, FATAL=0};
  template<class T> void write(log_level l, const T& v)
  { os << v; }
  void write(log_level l, const std::string& v);
  void write(log_level l, const char* v)
  { write(l, std::string(v)); }

//-----------------------------------------------------------------------
/// Underlying stream, can be null if no underlying stream.
//-----------------------------------------------------------------------

  virtual std::ostream* stream() = 0;

//-----------------------------------------------------------------------
/// Flush data to the log, at the given level.
//-----------------------------------------------------------------------

  virtual void flush(log_level l) = 0;

  void print(std::ostream& Os) {Os << "LogImp";}
protected:
  std::ostringstream os;
};

/****************************************************************//**
  This is a class that holds the level we are logging. It allows for
  things like Logger::debug() << "The answer is " << 5 << "\n".
*******************************************************************/

class LogHelper {
public:
  LogHelper(LogImp::log_level l, boost::shared_ptr<LogImp>& imp) 
  : l_(l), imp_(imp) {}
  template<class T> LogHelper& operator<<(T v) 
  {
    if(imp_) 
      imp_->write(l_, v);
    return *this;
  }
private:
  LogImp::log_level l_;
  boost::shared_ptr<LogImp> imp_;
};

/****************************************************************//**
  This is a simple logger. The logger depends on a specific
  implementation being set by set_implementation. If we don't have an
  implementation set, then the logger doesn't do anything.

  The logger is a singleton, there is just one global logger. You can
  directly access it through instance, but normally you just do things
  like "Logger::debug() << 'My debug message\n'". Data is flushed when
  a "\n" is encountered, you should end all messages with "\n".
*******************************************************************/

class Logger : public Printable<Logger> {
public:
//-----------------------------------------------------------------------
/// Set the implementation. It is perfectly legal for this to be a
/// null pointer, in that case we just don't send log messages
/// anywhere. 
//-----------------------------------------------------------------------

  static void set_implementation(const boost::shared_ptr<LogImp>& imp)
  { instance().imp_ = imp; }

//-----------------------------------------------------------------------
/// Set the implementation. It is perfectly legal for this to be a
/// null pointer, in that case we just don't send log messages
/// anywhere. 
//-----------------------------------------------------------------------

  static void set_implementation(LogImp* imp)
  { instance().imp_.reset(imp); }

//-----------------------------------------------------------------------
/// Underlying stream, can be null if no underlying stream.
//-----------------------------------------------------------------------

  static std::ostream* stream() 
  {
    if(instance().imp_)
      return instance().imp_->stream();
    else
      return 0;
  }
						    
//-----------------------------------------------------------------------
/// The logger instance.
//-----------------------------------------------------------------------

  static Logger& instance();
  static LogHelper log(LogImp::log_level l) 
  {return LogHelper(l, instance().imp_);}
  static LogHelper debug() {return log(LogImp::DEBUG);}
  static LogHelper info() {return log(LogImp::INFO);}
  static LogHelper warning() {return log(LogImp::WARNING);}
  static LogHelper error() {return log(LogImp::ERROR);}
  static LogHelper fatal() {return log(LogImp::FATAL);}
  void print(std::ostream& Os) {Os << "Logger";}
private:
  Logger() {}
  Logger(const Logger&) {}
  boost::shared_ptr<LogImp> imp_;
};

}
#endif
