#include "relative_humidity.h"
#include "ostream_pad.h"
#include "default_constant.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(RelativeHumidity)
.def(luabind::constructor<const boost::shared_ptr<Absorber>&, 
     const boost::shared_ptr<Temperature>&,
     const boost::shared_ptr<Pressure>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

RelativeHumidity::RelativeHumidity
(const boost::shared_ptr<Absorber>& Abs, 
 const boost::shared_ptr<Temperature>& Temp,
 const boost::shared_ptr<Pressure>& Press)
  : absorber(Abs), temp(Temp), press(Press)
{ 
  DefaultConstant cs;
  c = (cs.molar_weight_dry_air() / cs.molar_weight_water()).value;
}

boost::shared_ptr<RelativeHumidity> RelativeHumidity::clone() const
{
  return boost::shared_ptr<RelativeHumidity>(new RelativeHumidity(absorber, temp, press));
}


boost::shared_ptr<RelativeHumidity> 
RelativeHumidity::clone(const boost::shared_ptr<Absorber>& Abs, 
			const boost::shared_ptr<Temperature>& Temp,
			const boost::shared_ptr<Pressure>& Press) const
{
  return boost::shared_ptr<RelativeHumidity>(new RelativeHumidity(Abs, Temp, Press));
}

//-----------------------------------------------------------------------
/// Print to a string.
//-----------------------------------------------------------------------
void RelativeHumidity::print(std::ostream& Os)
{
  OstreamPad opad(Os, "    ");
  Os << "RelativeHumidity:\n";
  Os << "  Absorber:\n";
  opad << *absorber << "\n";
  opad.strict_sync();
  Os << "  Temperature:\n";
  opad << *temp << "\n";
  opad.strict_sync();
  Os << "  Pressure:\n";
  opad << *press << "\n";
  opad.strict_sync();
}

//-----------------------------------------------------------------------
/// Calculate specific humidity.
//-----------------------------------------------------------------------

ArrayAd<double, 1> RelativeHumidity::specific_humidity_grid() const
{
  blitz::Array<AutoDerivative<double>, 1> vgrid =
    absorber->absorber_vmr("H2O")->vmr_grid(*press).to_array();
  blitz::Array<AutoDerivative<double>, 1> shgrid(vgrid.rows());
  shgrid = vgrid / (c + vgrid);
  return ArrayAd<double, 1>(shgrid);
}

//-----------------------------------------------------------------------
/// Calculate relative humidity.
//-----------------------------------------------------------------------

ArrayAd<double, 1> RelativeHumidity::relative_humidity_grid() const
{
  // This comes Suniti. Not sure of a reference for this.
  blitz::Array<AutoDerivative<double>, 1> pgrid = 
    press->pressure_grid().convert(Unit("Pa")).value.to_array();
  blitz::Array<AutoDerivative<double>, 1> shgrid =
    specific_humidity_grid().to_array();
  blitz::Array<AutoDerivative<double>, 1> tgrid =
    temp->temperature_grid(*press).convert(Unit("K")).value.to_array();
  blitz::Array<AutoDerivative<double>, 1> res(pgrid.rows());
  res = 0.263 * pgrid * shgrid * exp(-17.67*(tgrid-273.16)/(tgrid-29.65));
  return ArrayAd<double, 1>(res);
}

//-----------------------------------------------------------------------
/// Relative humidity for each layer. This is just the average of the 
/// 2 levels.
//-----------------------------------------------------------------------

ArrayAd<double, 1> RelativeHumidity::relative_humidity_layer() const
{
  ArrayAd<double, 1> rh_lev = relative_humidity_grid();
  ArrayAd<double, 1> res(rh_lev.rows() - 1, rh_lev.number_variable());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = (rh_lev(i) + rh_lev(i+1)) / 2.0;
  return res;
}
