#include "temperature_met.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(TemperatureMet, Temperature)
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
			  const boost::shared_ptr<Pressure>&,
			  double, bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create an Temperature. 
//-----------------------------------------------------------------------

TemperatureMet::TemperatureMet
(const boost::shared_ptr<Meteorology>& Met_file,
 const boost::shared_ptr<Pressure>& Press,
 double Temp_offset,
 bool Temp_flag)
: TemperatureOffset(Press, Temp_offset, Temp_flag), met(Met_file)
{
}

// See base class for description of this function
boost::shared_ptr<Temperature> 
TemperatureMet::clone(const boost::shared_ptr<Pressure>& Press) const
{
  boost::shared_ptr<Temperature> res
    (new TemperatureMet(met, Press, coefficient()(0).value(),
			  used_flag_value()(0)));
  return res;
}

void TemperatureMet::print(std::ostream& Os) const 
{
  OstreamPad opad(Os, "    ");
  Os << "TemperatureMet:\n"
     << "  Temperature offset: " << temperature_offset() << "\n"
     << "  Retrieval flag:     " << (used_flag_value()(0) ? 
					"True\n" : "False\n")
     << "  Meteorology:\n";
  opad << *met << "\n";
  opad.strict_sync();
}
