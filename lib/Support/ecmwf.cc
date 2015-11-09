#include "ecmwf.h"
#include "log_interpolate.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
typedef Array<double, 1> (Ecmwf::*f1)(const Array<double, 1>& Pressure_level) const;
REGISTER_LUA_CLASS(Ecmwf)
.def("surface_pressure", &Ecmwf::surface_pressure)
.def("windspeed", &Ecmwf::windspeed)
.def("windspeed", &Ecmwf::windspeed_u)
.def("windspeed", &Ecmwf::windspeed_v)
.def("temperature", ((f1) &Ecmwf::temperature))
.def("h2o_vmr", ((f1) &Ecmwf::h2o_vmr))
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Return the H20 VMR. This is the specific_humidity converted to a
/// volume mixing ratio.
//-----------------------------------------------------------------------

blitz::Array<double, 1> Ecmwf::h2o_vmr(const blitz::Array<double, 1>& Pressure_level) const
{
  Array<double, 1> s = specific_humidity(Pressure_level);
  Array<double, 1> vmr(s.shape());
  vmr = s / (1 - s) * OldConstant::molar_weight_dry_air / 
    OldConstant::molar_weight_water;
  return vmr;
}

ArrayAd<double, 1> Ecmwf::h2o_vmr(const ArrayAd<double, 1>& Pressure_level) const
{
  Array<AutoDerivative<double>, 1> s = 
    specific_humidity(Pressure_level).to_array();
  Array<AutoDerivative<double>, 1> vmr(s.shape());
  vmr = s / (1 - s) * OldConstant::molar_weight_dry_air / 
    OldConstant::molar_weight_water;
  return vmr;
}

//-----------------------------------------------------------------------
/// The temperature and specific_humidity reading is identical except
/// for the field read. So we have a generic function here with the
/// field name passed in.
//-----------------------------------------------------------------------

blitz::Array<double, 1> Ecmwf::read_and_interpolate(const std::string& Field,
	    const blitz::Array<double, 1>& Pressure_level) const
{
  Array<double, 1> p, t;
  read(Field, p, t);
  LogLogInterpolate<double, double> tint(p.begin(), p.end(), t.begin());
  Array<double, 1> tres(Pressure_level.shape());
  for(int i = 0; i < Pressure_level.rows(); ++i)
    tres(i) = tint(Pressure_level(i));
  return tres;
}

ArrayAd<double, 1> Ecmwf::read_and_interpolate(const std::string& Field,
	    const ArrayAd<double, 1>& Pressure_level) const
{
  firstIndex i1; secondIndex i2;
  Array<double, 1> p, t;
  read(Field, p, t);
  std::vector<AutoDerivative<double> > p2, t2;
  for(int i = 0; i < p.rows(); ++i) {
    p2.push_back(p(i));
    t2.push_back(t(i));
  }
  LogLogInterpolate<AutoDerivative<double>, AutoDerivative<double> > 
    tint(p2.begin(), p2.end(), t2.begin());
  Array<AutoDerivative<double>, 1> tres(Pressure_level.rows());
  for(int i = 0; i < Pressure_level.rows(); ++i)
    tres(i) = tint(Pressure_level(i));
  return tres;
}
