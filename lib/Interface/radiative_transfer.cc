#include "radiative_transfer.h"
#include "logger.h"
using namespace FullPhysics;

AccumulatedTimer RadiativeTransfer::timer("RT");

using namespace FullPhysics;
#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(RadiativeTransfer)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Helper routine, creates a progress meter. This will return 0 if we
/// aren't logging, or if we don't have enough points to bother with.
//-----------------------------------------------------------------------
boost::shared_ptr<boost::progress_display> 
RadiativeTransfer::progress_display(const blitz::Array<double, 1>& wn) const
{
  boost::shared_ptr<boost::progress_display> res;
  if(wn.size() > 100 &&
     Logger::stream()) {
    Logger::info() << "RT Progress\n";
    res.reset(new boost::progress_display(wn.size(), *Logger::stream()));
  }
  return res;
}
