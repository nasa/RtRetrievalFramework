#include "ground_coxmunk_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<RegisterOutputBase> ground_cm_output_create(boost::shared_ptr<Ground>& ground)
{
    return boost::shared_ptr<RegisterOutputBase>
        (new GroundCoxmunkOutput(boost::dynamic_pointer_cast<GroundCoxmunk>(ground)));
}

REGISTER_LUA_DERIVED_CLASS(GroundCoxmunkOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<GroundCoxmunk>&>())
.scope
[
    luabind::def("create", &ground_cm_output_create)
]
REGISTER_LUA_END()
#endif

double windspeed(boost::shared_ptr<GroundCoxmunk>& Coxmunk)
{
  return Coxmunk->coefficient()(0).value();
}

double windspeed_uncert(boost::shared_ptr<GroundCoxmunk>& Coxmunk)
{
  Array<double, 2> cov = Coxmunk->statevector_covariance();
  if(cov.rows() > 0 and cov.rows() > 0) {
    return (cov(0, 0) < 0 ? 0.0 : sqrt(cov(0, 0)));
  } else {
    return 0.0;
  }
}

void GroundCoxmunkOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  boost::shared_ptr<GroundCoxmunk> cm_freeze = 
    boost::dynamic_pointer_cast<GroundCoxmunk>(coxmunk->clone());

  { boost::function<double ()> f = boost::bind(&windspeed, cm_freeze);
    out->register_data_source("/RetrievalResults/wind_speed_apriori", f); }

}

void GroundCoxmunkOutput::register_output(const boost::shared_ptr<Output>& out) const
{

  { boost::function<double ()> f = boost::bind(&windspeed, coxmunk);
    out->register_data_source("/RetrievalResults/wind_speed", f); }

  { boost::function<double ()> f = boost::bind(&windspeed_uncert, coxmunk);
    out->register_data_source("/RetrievalResults/wind_speed_uncertainty", f); }

  out->register_data_source("/RetrievalResults/surface_type", surface_type.c_str()); 
}
