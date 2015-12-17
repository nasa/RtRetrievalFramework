// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "fp_logger.h"
%}
%base_import(logger)
%fp_shared_ptr(FullPhysics::FpLogger);

namespace FullPhysics {
class FpLogger : public LogImp {
public:
  FpLogger(int Verbosity_level = LogImp::INFO);
  virtual void flush(log_level l);
  %extend {
    static void turn_on_logger()
    {
      FullPhysics::Logger::set_implementation(new FullPhysics::FpLogger);
    }
    static void turn_off_logger()
    {
      FullPhysics::Logger::set_implementation(0);
    }
  }
};
}

%init %{
  FullPhysics::Logger::set_implementation(new FullPhysics::FpLogger);
%}

