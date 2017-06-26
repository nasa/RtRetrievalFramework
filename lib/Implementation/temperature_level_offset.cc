#include "temperature_level_offset.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(TemperatureLevelOffset, Temperature)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
			  const blitz::Array<double, 1>&,
			  double, bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create an Temperature. 
//-----------------------------------------------------------------------

TemperatureLevelOffset::TemperatureLevelOffset
(const boost::shared_ptr<Pressure>& Press,
 const blitz::Array<double, 1>& Temp_levels,
 double Temp_offset,
 bool Temp_flag)
: TemperatureOffset(Press, Temp_offset, Temp_flag), temp_levels(Temp_levels)
{
}

// See base class for description of this function

boost::shared_ptr<Temperature> 
TemperatureLevelOffset::clone(const boost::shared_ptr<Pressure>& Press) const
{
  boost::shared_ptr<Temperature> res
    (new TemperatureLevelOffset(Press, temp_levels, coefficient()(0).value(),
			  used_flag_value()(0)));
  return res;
}

void TemperatureLevelOffset::print(std::ostream& Os) const 
{
  OstreamPad opad(Os, "    ");
  Os << "TemperatureLevelOffset:\n"
     << "  Temperature offset: " << coefficient()(0) << "\n"
     << "  Retrieval flag:     " << (used_flag_value()(0) ? 
					"True\n" : "False\n")
     << "  Level Initial Guess:\n";
  opad << temp_levels << "\n";
  opad.strict_sync();
}
