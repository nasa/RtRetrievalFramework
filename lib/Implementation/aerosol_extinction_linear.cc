#include "aerosol_extinction_linear.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolExtinctionLinear, AerosolExtinction)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
			  const blitz::Array<bool, 1>&, 
			  const blitz::Array<double, 1>&,
			  const std::string&>())
REGISTER_LUA_END()
#endif

// See base class for description
boost::shared_ptr<AerosolExtinction> AerosolExtinctionLinear::clone
(const boost::shared_ptr<Pressure>& Pres) const
{
  return boost::shared_ptr<AerosolExtinction>
    (new AerosolExtinctionLinear(Pres, used_flag, coeff.value(), 
				 aerosol_name()));
}

void AerosolExtinctionLinear::calc_aerosol_extinction() const
{
  aext.resize(coeff.rows(), coeff.number_variable());
  for(int i = 0; i < coeff.rows(); ++i)
    aext(i) = coeff(i);
}

void AerosolExtinctionLinear::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "    ");
  Os << "AerosolExtinctionLinear:\n"
     << "  Coefficient:\n";
  opad << coeff.value() << "\n";
  opad.strict_sync();
  Os << "  Retrieval flag:\n";
  opad << used_flag << "\n";
  opad.strict_sync();
}
