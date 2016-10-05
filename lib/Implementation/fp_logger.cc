#include "fp_logger.h"
#include "fp_exception.h"
#include <iostream>

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(FpLogger, LogImp)
.def(luabind::constructor<int>())
REGISTER_LUA_END()
#endif


//-----------------------------------------------------------------------
/// Flush data to log
//-----------------------------------------------------------------------

void FpLogger::flush(log_level l)
{
  std::string s = os.str();
  os.str("");
  if((int) l > verbosity_level_)
    return;
  switch(l) {
  case LogImp::DEBUG:
    *stream() << "DEBUG:   " << s;
    stream()->flush();
    break;
  case LogImp::INFO:
    *stream() << "INFO:    " << s;
    stream()->flush();
    break;
  case LogImp::WARNING:
    *stream() << "WARNING: " << s;
    stream()->flush();
    break;
  case LogImp::ERROR:
    std::cerr << "ERROR:   " << s;
    break;
  case LogImp::FATAL:
    std::cerr << "FATAL:   " << s;
    break;
  default:
    throw Exception("Unknown log level");
  }
}
