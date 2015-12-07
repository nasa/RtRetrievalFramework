#include "relative_humidity.h"
#include "ostream_pad.h"
#include "default_constant.h"
using namespace FullPhysics;

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
/// Calculate relative humidity.
//-----------------------------------------------------------------------

ArrayAd<double, 1> RelativeHumidity::relative_humidity_grid() const
{
  // This comes Suniti. Not sure of a reference for this.
  blitz::Array<AutoDerivative<double>, 1> pgrid = 
    press->pressure_grid().convert(Unit("Pa")).value.to_array();
  blitz::Array<AutoDerivative<double>, 1> vgrid =
    absorber->absorber_vmr("H2O")->vmr_grid(*press).to_array();
  blitz::Array<AutoDerivative<double>, 1> shgrid(vgrid.rows());
  shgrid = vgrid / (c + vgrid);
  blitz::Array<AutoDerivative<double>, 1> tgrid =
    temp->temperature_grid(*press).convert(Unit("K")).value.to_array();
  blitz::Array<AutoDerivative<double>, 1> res(pgrid.rows());
  res = 0.263 * pgrid * shgrid * exp(-17.67*(tgrid-273.16)/(tgrid-29.65));
  return ArrayAd<double, 1>(res);
}
