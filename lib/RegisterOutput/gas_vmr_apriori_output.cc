#include "gas_vmr_apriori_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(GasVmrAprioriOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<GasVmrApriori>&>())
REGISTER_LUA_END()
#endif

double tropopause_altitude(const boost::shared_ptr<GasVmrApriori>& ga)
{
  return ga->tropopause_altitude().convert(units::m).value;
}

double tropopause_pressure(const boost::shared_ptr<GasVmrApriori>& ga)
{
  return ga->tropopause_pressure();
}

void GasVmrAprioriOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
}

void GasVmrAprioriOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    { boost::function<double ()> f = boost::bind(&tropopause_altitude, gas_apriori); 
      out->register_data_source("/RetrievalResults/tropopause_altitude", f); }

    { boost::function<double ()> f = boost::bind(&tropopause_pressure, gas_apriori); 
      out->register_data_source("/RetrievalResults/tropopause_pressure", f); }
}
