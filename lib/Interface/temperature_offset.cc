#include "temperature_offset.h"
#include "linear_interpolate.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(TemperatureOffset, Temperature)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create an Temperature Offset. 
//-----------------------------------------------------------------------

TemperatureOffset::TemperatureOffset(const boost::shared_ptr<Pressure>& Press,
				     double Temp_offset,
				     bool Temp_flag)
{
  Array<bool, 1> flag(1);
  Array<double, 1> val(flag.shape());
  flag(0) = Temp_flag;
  val(0) = Temp_offset;
  init(val, flag, Press, false);
}

//-----------------------------------------------------------------------
/// This calculates temperature grid to use for layer retrieval. 
//-----------------------------------------------------------------------

void TemperatureOffset::calc_temperature_grid() const
{
  blitz::Array<double, 1> temp_profile( temperature_profile() );
  blitz::Array<double, 1> press_profile( pressure_profile() );

  std::vector<AutoDerivative<double> > plist;
  std::vector<AutoDerivative<double> > tlist;
  if (press_profile.rows() != temp_profile.rows()) {
    std::stringstream err_msg;
    err_msg << "Size of pressure grid: "
	    << press_profile.rows()
	    << " != size of temperature levels: "
	    << temp_profile.rows();
    throw Exception(err_msg.str());
  }
  for(int i = 0; i < press_profile.rows(); ++i) {
    AutoDerivative<double> t2 = temp_profile(i) + coefficient()(0);
    tlist.push_back(t2);
    plist.push_back(press_profile(i));
  }
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  boost::shared_ptr<lin_type> lin
    (new lin_type(plist.begin(), plist.end(), tlist.begin()));
  tgrid = boost::bind(&lin_type::operator(), lin, _1);
}

// See base class for description of this function.
std::string TemperatureOffset::state_vector_name_i(int i) const
{
  return "Temperature Offset (Kelvin)";
}

//-----------------------------------------------------------------------
/// Uncertainty of temperature offset.
//-----------------------------------------------------------------------

double TemperatureOffset::temperature_offset_uncertainty() const
{
  if(!used_flag_value()(0) ||
     cov.rows() == 0 ||
     cov(0,0) < 0)
    return 0.0;
  return sqrt(cov(0,0));
}

void TemperatureOffset::print(std::ostream& Os) const 
{
  Os << "TemperatureOffset\n";
}
