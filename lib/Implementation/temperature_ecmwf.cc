#include "temperature_ecmwf.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(TemperatureEcmwf, Temperature)
.def(luabind::constructor<const boost::shared_ptr<Ecmwf>&,
			  const boost::shared_ptr<Pressure>&,
			  double, bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create an Temperature. 
//-----------------------------------------------------------------------

TemperatureEcmwf::TemperatureEcmwf
(const boost::shared_ptr<Ecmwf>& Ecmwf_file,
 const boost::shared_ptr<Pressure>& Press,
 double Temp_offset,
 bool Temp_flag)
: TemperatureOffset(Press, Temp_offset, Temp_flag), ecmwf(Ecmwf_file)
{
}

// See base class for description of this function
boost::shared_ptr<Temperature> 
TemperatureEcmwf::clone(const boost::shared_ptr<Pressure>& Press) const
{
  boost::shared_ptr<Temperature> res
    (new TemperatureEcmwf(ecmwf, Press, coefficient()(0).value(),
			  used_flag_value()(0)));
  return res;
}

void TemperatureEcmwf::print(std::ostream& Os) const 
{
  OstreamPad opad(Os, "    ");
  Os << "TemperatureEcmwf:\n"
     << "  Temperature offset: " << temperature_offset() << "\n"
     << "  Retrieval flag:     " << (used_flag_value()(0) ? 
					"True\n" : "False\n")
     << "  ECMWF:\n";
  opad << *ecmwf << "\n";
  opad.strict_sync();
}
