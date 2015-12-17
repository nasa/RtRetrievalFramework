#include "absorber_vmr_fixed_level_scaled.h"
#include "ostream_pad.h"
#include "linear_interpolate.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrFixedLevelScaled, AbsorberVmr)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
			  const boost::shared_ptr<PressureLevelInput>&,
			  const blitz::Array<double, 1>&,
			  bool,
			  double,
			  const std::string&>())
REGISTER_LUA_END()
#endif

AbsorberVmrFixedLevelScaled:: AbsorberVmrFixedLevelScaled
(const boost::shared_ptr<Pressure>& Press,
 const boost::shared_ptr<PressureLevelInput>& Press_level,
 const blitz::Array<double, 1>& Vmr,
 bool Used_flag,
 double Scale,			      
 const std::string& Gas_name)
: press_level(Press_level),
  vmr0(Vmr.copy())
{
  Array<bool, 1> uflag(1);
  Array<double, 1> val(1);
  val(0) = Scale;
  uflag(0) = Used_flag;
  init(Gas_name, val, uflag, Press, false);
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrFixedLevelScaled::clone() const
{
  boost::shared_ptr<Pressure> pressure_clone = press->clone();
  return clone(pressure_clone);
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrFixedLevelScaled::clone
(const boost::shared_ptr<Pressure>& Press) const
{
  return boost::shared_ptr<AbsorberVmr>
    (new AbsorberVmrFixedLevelScaled(Press, press_level, vmr0, 
				     used_flag(0), scale_factor(), gas_name()));
}

void AbsorberVmrFixedLevelScaled::calc_vmr() const
{
  AutoDerivative<double> p = log(press->surface_pressure().value);
  std::vector<AutoDerivative<double> > plist;
  std::vector<AutoDerivative<double> > vmrlist;
  for(int i = 0; i < press->pressure_grid().rows() - 1; ++i) {
    vmrlist.push_back(vmr0(i) * coeff(0));
    plist.push_back(press->pressure_grid()(i).value);
  }
  int i = press->pressure_grid().rows() - 1;
  double p1 = log(press_level->pressure_level()(i - 1));
  double p2 = log(press_level->pressure_level()(i));
  AutoDerivative<double> v = coeff(0) *
    (vmr0(i - 1) + (p - p1) * (vmr0(i) - vmr0(i - 1)) / (p2 - p1));
  plist.push_back(press->surface_pressure().value);
  vmrlist.push_back(v);
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  boost::shared_ptr<lin_type> lin
    (new lin_type(plist.begin(), plist.end(), vmrlist.begin()));
  vmr = boost::bind(&lin_type::operator(), lin, _1);
}

void AbsorberVmrFixedLevelScaled::print(std::ostream& Os) const
{ 
  OstreamPad opad(Os, "    ");
  Os << "AbsorberVmrFixedLevelScaled:\n"
     << "  Gas name: " << gas_name() << "\n"
     << "  Scale:    " << scale_factor() << "\n"
     << "  VMR:\n";
  opad << vmr0 << "\n";
  opad.strict_sync();
}
