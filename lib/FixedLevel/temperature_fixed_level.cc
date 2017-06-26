#include "temperature_fixed_level.h"
#include "fp_exception.h"
#include "pressure_fixed_level.h"
#include "ostream_pad.h"
#include "linear_interpolate.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(TemperatureFixedLevel, Temperature)
.def(luabind::constructor<const blitz::Array<bool, 1>&, bool, 
			  const blitz::Array<double, 1>&,
			  double, 
			  const boost::shared_ptr<Pressure>&,
			  const boost::shared_ptr<PressureLevelInput>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create an Temperature. 
//-----------------------------------------------------------------------

TemperatureFixedLevel::TemperatureFixedLevel
(const blitz::Array<bool, 1>& Flag_temp, 
 bool Flag_offset, const blitz::Array<double, 1>& Temp,
 double T_offset, 
 const boost::shared_ptr<Pressure>& Press,
 const boost::shared_ptr<PressureLevelInput>& Press_level)
: press_level(Press_level)
{
  if(Flag_temp.rows() != Temp.rows())
    throw Exception("Flag_temp and Temp need to be the same size");
  Range trange(1, Flag_temp.rows());
  Array<bool, 1> flag(Flag_temp.rows() + 1);
  Array<double, 1> val(flag.shape());
  flag(0) = Flag_offset;
  val(0) = T_offset;
  flag(trange) = Flag_temp;
  val(trange) = Temp;
  // When marking, only include the temperature values if they are below the
  // number of levels given by Press.
  init(val, flag, Press, true, 1);
}

//-----------------------------------------------------------------------
/// This calculates temperature grid to use for layer retrieval. This
/// has the same number of layers are press->number_layer(), and we
/// interpolate the temperature() at the fixed levels to the surface
/// pressure. 
///
/// As the surface pressure is changed, the size of the pressure
/// grid can change (e.g., surface pressure crosses one of the
/// pressure levels).
//-----------------------------------------------------------------------

void TemperatureFixedLevel::calc_temperature_grid() const
{
  AutoDerivative<double> p = log(press->surface_pressure().value);
  std::vector<AutoDerivative<double> > plist;
  std::vector<AutoDerivative<double> > tlist;
  for(int i = 0; i < press->pressure_grid().rows() - 1; ++i) {
    plist.push_back(press->pressure_grid()(i).value);
    tlist.push_back(temperature_levels()(i));
  }
  int i = press->pressure_grid().rows() - 1;
  double p1 = log(press_level->pressure_level()(i - 1));
  AutoDerivative<double> t1, t2;
  t1 = temperature_levels()(i - 1);
  double p2 = log(press_level->pressure_level()(i));
  t2 = temperature_levels()(i);
  AutoDerivative<double> t = t1 + (p - p1) * (t2 - t1) / (p2 - p1);
  plist.push_back(press->surface_pressure().value);
  tlist.push_back(t);
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  boost::shared_ptr<lin_type> lin
    (new lin_type(plist.begin(), plist.end(), tlist.begin()));
  tgrid = boost::bind(&lin_type::operator(), lin, _1);
}

// See base class for description of this function

boost::shared_ptr<Temperature> 
TemperatureFixedLevel::clone(const boost::shared_ptr<Pressure>& Press) const
{
  boost::shared_ptr<Temperature> res
    (new TemperatureFixedLevel(used_flag_value()(temperature_range()), 
			       used_flag_value()(0),
			       coefficient()(temperature_range()).value(),
			       coefficient()(0).value(), Press, press_level));
  return res;
}

// See base class for description of this function.

std::string TemperatureFixedLevel::state_vector_name_i(int i) const
{
  if(i ==0)
    return "Temperature Offset (Kelvin)";
  return "Temperature (Kelvin) for Pressure Level " + 
    boost::lexical_cast<std::string>(i);
}

//-----------------------------------------------------------------------
/// Return the temperature on the fixed levels (which may include values
/// from below the surface). This is in Kelvin.
//-----------------------------------------------------------------------

ArrayAd<double, 1> 
TemperatureFixedLevel::temperature_levels() const
{ 
  ArrayAd<double, 1> 
    res(coefficient()(temperature_range()).copy());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = res(i) + coefficient()(0);
  return res;
}

//-----------------------------------------------------------------------
/// Uncertainty of temperature offset.
//-----------------------------------------------------------------------

double TemperatureFixedLevel::temperature_offset_uncertainty() const
{
  if(!used_flag_value()(0) ||
     cov.rows() == 0 ||
     cov(0,0) < 0)
    return 0.0;
  return sqrt(cov(0,0));
}

void TemperatureFixedLevel::print(std::ostream& Os) const 
{
  OstreamPad opad(Os, "    ");
  Os << "TemperatureFixedLevel:\n"
     << "  Temperature offset:    " << temperature_offset() << "\n"
     << "  Retrieval flag offset: " << (used_flag_value()(0) ? 
					"True\n" : "False\n")
     << "  Temperature:\n";
  ArrayAd<double, 1> 
    temp(coefficient()(temperature_range()));
  opad << temp.value() << "\n";
  opad.strict_sync();
  Os << "  Temperature flag:\n";
  opad << used_flag_value()(temperature_range()) << "\n";
  opad.strict_sync();
}
