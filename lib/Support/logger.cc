#include "logger.h"
#include "fp_exception.h"

using namespace FullPhysics;

using namespace FullPhysics;
#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(LogImp)
.enum_("log_level")
[
 luabind::value("DEBUG", LogImp::DEBUG),
 luabind::value("INFO", LogImp::INFO),
 luabind::value("WARNING", LogImp::WARNING),
 luabind::value("ERROR", LogImp::ERROR),
 luabind::value("FATAL", LogImp::FATAL)
]
REGISTER_LUA_END()


typedef void (*f1)(const boost::shared_ptr<LogImp>& imp);
REGISTER_LUA_CLASS(Logger)
.scope
[
 luabind::def("set_implementation", ((f1) &Logger::set_implementation))
]
REGISTER_LUA_END()
#endif


//-----------------------------------------------------------------------
/// Return the instance of the Logger.
//-----------------------------------------------------------------------

Logger& Logger::instance()
{
  static Logger l;
  return l;
}

//-----------------------------------------------------------------------
/// Write to a log.
//-----------------------------------------------------------------------

void LogImp::write(log_level l, const std::string& v)
{
  os << v;
  if(v.find("\n") != std::string::npos)
    flush(l);
}

extern "C" {
  void lg_write_log_wrap(const int* level, const char* s, const int* s_len);
}

void lg_write_log_wrap(const int* level, const char* s, const int* s_len)
{
  Logger::log(LogImp::log_level(*level)) << std::string(s,*s_len) << "\n";

  // Throw exception here since a FATAL from fortran code really indicates
  // that codes should have ended
  if (*level == LogImp::FATAL)
    throw Exception(std::string(s,*s_len));
}
