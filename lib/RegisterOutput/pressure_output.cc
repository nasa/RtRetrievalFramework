#include "pressure_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(PressureOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
     const boost::shared_ptr<StateVector>& >())
REGISTER_LUA_END()
#endif

class PressureOutputHelper {
public:
  PressureOutputHelper(const boost::shared_ptr<Pressure>& P,
		       const boost::shared_ptr<StateVector>& Sv)
    : p(P), sv(Sv) {}
  double surface_pressure_uncertainty() const
  {
    // There is no covariance value loaded so we can not return an uncertainty value
    // Must be using forward model mode
    if (sv->state_covariance().rows() == 0 or sv->state_covariance().cols() == 0)
        return 0;

    firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
    Array<double, 1> dsurf_dstate = p->surface_pressure().value.gradient();
    Array<double, 1> t(dsurf_dstate.rows());
    t = sum(sv->state_covariance()(i1, i2) * dsurf_dstate(i2), i2);
    double t2 = sum(dsurf_dstate * t);
    return (t2 > 0 ? sqrt(t2) : 0);
  }
  Array<double, 1> pressure_grid_value() const
  {
    return p->pressure_grid().convert(units::Pa).value.value();
  }
private:
  boost::shared_ptr<Pressure> p;
  boost::shared_ptr<StateVector> sv;
};

// See base class for description

void PressureOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the pressure state
  boost::shared_ptr<Pressure> pfreeze = p->clone();
  out->register_data_source("/RetrievalResults/surface_pressure_apriori_fph", 
			   &Pressure::surface_pressure_value, pfreeze);

  boost::shared_ptr<StateVector> sv_zero;
  boost::shared_ptr<PressureOutputHelper> 
    phelp(new PressureOutputHelper(pfreeze, sv_zero));
  out->register_data_source("/RetrievalResults/vector_pressure_levels_apriori",
			   &PressureOutputHelper::pressure_grid_value, phelp);
}

void PressureOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source("/RetrievalResults/surface_pressure_fph", 
			   &Pressure::surface_pressure_value, p);
  boost::shared_ptr<PressureOutputHelper> 
    phelp(new PressureOutputHelper(p, sv));
  out->register_data_source("/RetrievalResults/surface_pressure_uncert_fph", 
		   &PressureOutputHelper::surface_pressure_uncertainty, phelp);
  out->register_data_source("/RetrievalResults/vector_pressure_levels",
			   &PressureOutputHelper::pressure_grid_value, phelp);
  out->register_data_source("/RetrievalResults/num_active_levels",
			   &Pressure::number_level, p);
}

