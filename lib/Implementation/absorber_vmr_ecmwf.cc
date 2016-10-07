#include "absorber_vmr_ecmwf.h"
#include "ostream_pad.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrEcmwf, AbsorberVmr)
.def(luabind::constructor<const boost::shared_ptr<Ecmwf>&,
			  const boost::shared_ptr<Pressure>&,
			  double, 
			  bool,
			  const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

AbsorberVmrEcmwf::AbsorberVmrEcmwf
(const boost::shared_ptr<Ecmwf>& Ecmwf_file,
 const boost::shared_ptr<Pressure>& Press,
 double Scale,                         
 bool Scale_flag,
 const std::string& Gas_name)
: AbsorberVmrScaled(Press, Scale, Scale_flag, Gas_name), ecmwf(Ecmwf_file)
{
  std::string gname = gas_name();
  boost::to_upper(gname);
  if (gname != "H2O" &&
      gname != "O3") {
    std::stringstream err_msg; 
    err_msg << "Only H2O and O3 is supported by AbsorberVmrEcmwf, unknown absorber: " 
	    << gas_name();
    throw Exception(err_msg.str());
  }
}

blitz::Array<double, 1> AbsorberVmrEcmwf::specific_humidity_ecmwf() const
{ 
  blitz::Array<double, 1> s, p;
  ecmwf->specific_humidity_grid(p, s);
  return s;
}

blitz::Array<double, 1> AbsorberVmrEcmwf::vmr_profile() const
{ 
  if(gas_name() == "H2O") {
    blitz::Array<double, 1> s( specific_humidity_ecmwf() );
    blitz::Array<double, 1> res( s / (1 - s) * OldConstant::molar_weight_dry_air / 
				 OldConstant::molar_weight_water );
    return res;
  }
  if(gas_name() == "O3") {
    blitz::Array<double, 1> s, p;
    ecmwf->ozone_mmr_grid(p, s);
    blitz::Array<double, 1> res(s * OldConstant::molar_weight_dry_air / 
				 OldConstant::molar_weight_ozone);
    return res;
  }
}

blitz::Array<double, 1> AbsorberVmrEcmwf::pressure_profile() const
{ 
  blitz::Array<double, 1> s, p;
  ecmwf->specific_humidity_grid(p, s);
  return p;
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrEcmwf::clone
(const boost::shared_ptr<Pressure>& Press) const
{
  return boost::shared_ptr<AbsorberVmr>
    (new AbsorberVmrEcmwf(ecmwf, Press, coeff(0).value(),used_flag(0),
			  gas_name()));
}

void AbsorberVmrEcmwf::print(std::ostream& Os) const
{ 
  OstreamPad opad(Os, "    ");
  Os << "AbsorberVmrEcmwf:\n"
     << "  Gas name:       " << gas_name() << "\n"
     << "  Scale:          " << scale_factor() << "\n"
     << "  Retrieval flag: " << (used_flag_value()(0) ? 
					"True\n" : "False\n")
     << "  ECMWF:\n";
  opad << *ecmwf << "\n";
  opad.strict_sync();
}
