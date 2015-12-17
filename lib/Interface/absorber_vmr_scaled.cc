#include "absorber_vmr_scaled.h"
#include "linear_interpolate.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrScaled, AbsorberVmr)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

AbsorberVmrScaled::AbsorberVmrScaled
(const boost::shared_ptr<Pressure>& Press,
 double Scale,                         
 bool Scale_flag,
 const std::string& Gas_name)
{
  Array<bool, 1> flag(1);
  Array<double, 1> val(flag.shape());
  flag(0) = Scale_flag;
  val(0) = Scale;
  init(Gas_name, val, flag, Press, false);
}

void AbsorberVmrScaled::calc_vmr() const
{
  blitz::Array<double, 1> v_profile( vmr_profile() );
  blitz::Array<double, 1> press_profile( pressure_profile() );

  std::vector<AutoDerivative<double> > plist;
  std::vector<AutoDerivative<double> > vlist;
  if (press_profile.rows() != v_profile.rows()) {
    std::stringstream err_msg;
    err_msg << "Size of pressure grid: "
	    << press_profile.rows()
	    << " != size of vmr levels: "
	    << v_profile.rows();
    throw Exception(err_msg.str());
  }
  for(int i = 0; i < press_profile.rows(); ++i) {
    AutoDerivative<double> t2 = v_profile(i) * coefficient()(0);
    vlist.push_back(t2);
    plist.push_back(press_profile(i));
  }
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  boost::shared_ptr<lin_type> lin
    (new lin_type(plist.begin(), plist.end(), vlist.begin()));
  vmr = boost::bind(&lin_type::operator(), lin, _1);
}

void AbsorberVmrScaled::print(std::ostream& Os) const
{ 
  Os << "AbosorberVmrScaled\n";
}
