// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "logger.h"
%}

%fp_shared_ptr(FullPhysics::LogImp);

namespace FullPhysics {
class LogImp {
public:
  virtual ~LogImp();
  enum log_level {DEBUG = 4, INFO=3, WARNING=2, ERROR=1, FATAL=0};
  void write(log_level l, const std::string& v);
  virtual void flush(log_level l) = 0;
};

%nodefaultctor Logger;
class Logger {
public:
  static void set_implementation(const boost::shared_ptr<LogImp>& imp);
};
}

