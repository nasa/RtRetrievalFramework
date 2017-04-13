#include "pressure_fixed_level_output.h"
#include "fill_value.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> press_create
(const boost::shared_ptr<Pressure>& P, const boost::shared_ptr<StateVector>& Sv)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new PressureFixedLevelOutput
     (boost::dynamic_pointer_cast<PressureFixedLevel>(P), Sv));
}
REGISTER_LUA_DERIVED_CLASS(PressureFixedLevelOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &press_create)
]
REGISTER_LUA_END()
#endif

class PressureFixedLevelOutputHelper {
public:
  PressureFixedLevelOutputHelper(const boost::shared_ptr<Pressure>& P,
		       const boost::shared_ptr<StateVector>& Sv)
    : p(P), sv(Sv) {}
  double surface_pressure_uncertainty() const
  {
    firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
    Array<double, 1> dsurf_dstate = p->surface_pressure().value.gradient();
    Array<double, 1> t(dsurf_dstate.rows());
    t = sum(sv->state_covariance()(i1, i2) * dsurf_dstate(i2), i2);
    double t2 = sum(dsurf_dstate * t);
    return (t2 > 0 ? sqrt(t2) : 0);
  }
  Array<double, 1> pressure_grid_value() const
  {
    return p->pressure_grid().value.value();
  }
private:
  boost::shared_ptr<Pressure> p;
  boost::shared_ptr<StateVector> sv;
};

void PressureFixedLevelOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the pressure state
  boost::shared_ptr<Pressure> pfreeze = p->clone();
  out->register_data_source("/RetrievalResults/surface_pressure_apriori_fph", 
			   &Pressure::surface_pressure_value, pfreeze);
}

// See base class for description

void PressureFixedLevelOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source("/RetrievalResults/surface_pressure_fph", 
			   &Pressure::surface_pressure_value, 
			   boost::dynamic_pointer_cast<Pressure>(p));
  boost::shared_ptr<PressureFixedLevelOutputHelper> 
    phelp(new PressureFixedLevelOutputHelper(p, sv));
  out->register_data_source("/RetrievalResults/surface_pressure_uncert_fph", 
   &PressureFixedLevelOutputHelper::surface_pressure_uncertainty, phelp);
  out->register_data_source("/RetrievalResults/num_active_levels", 
			   &PressureFixedLevel::number_active_level, p);
  out->register_data_source_pad("/RetrievalResults/vector_pressure_levels", 
			       &PressureFixedLevel::pressure_active_levels, p,
			       p->max_number_level(), fill_value<double>());
}

