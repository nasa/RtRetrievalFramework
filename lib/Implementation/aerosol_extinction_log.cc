#include "aerosol_extinction_log.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolExtinctionLog, AerosolExtinction)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
			  const blitz::Array<bool, 1>&, 
			  const blitz::Array<double, 1>&,
			  const std::string&>())
REGISTER_LUA_END()
#endif

// See base class for description
boost::shared_ptr<AerosolExtinction> AerosolExtinctionLog::clone
(const boost::shared_ptr<Pressure>& Pres) const
{
  return boost::shared_ptr<AerosolExtinction>
    (new AerosolExtinctionLog(Pres, used_flag, coeff.value(), aerosol_name()));
}

void AerosolExtinctionLog::calc_aerosol_extinction() const
{
  aext.resize(coeff.rows(), coeff.number_variable());
  for(int i = 0; i < coeff.rows(); ++i)
    aext(i) = exp(coeff(i));
}

void AerosolExtinctionLog::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "    ");
  Os << "AerosolExtinctionLog:\n"
     << "  Coefficient:\n";
  opad << coeff.value() << "\n";
  opad.strict_sync();
  Os << "  Retrieval flag:\n";
  opad << used_flag << "\n";
  opad.strict_sync();
}
