#ifndef FP_LOGGER_H
#define FP_LOGGER_H
#include "logger.h"
#include <iostream>

namespace FullPhysics {
/****************************************************************//**
  This is the implementation of the Logger used for the Full Physics 
  program. This just writes to stdout or stderr, filtering by the
  level, and adding in a leading label (e.g., "INFO").
*******************************************************************/
class FpLogger : public LogImp {
public:
//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

  FpLogger(int Verbosity_level = LogImp::DEBUG)
    : verbosity_level_(Verbosity_level) { }
  virtual ~FpLogger() {}
  virtual void flush(log_level l);
  virtual std::ostream* stream() {return &std::cout;}
private:
  int verbosity_level_;
};
}

#endif
