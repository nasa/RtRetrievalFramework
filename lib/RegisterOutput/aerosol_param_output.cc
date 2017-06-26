#include <boost/algorithm/string.hpp>
#include "aerosol_param_output.h"
#include "aerosol_extinction.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> aco_create
(const boost::shared_ptr<AerosolExtinction>& A)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new AerosolParamOutput
     (boost::dynamic_pointer_cast<AerosolExtinctionImpBase>(A)));
}
boost::shared_ptr<RegisterOutputBase> aco_create2
(const boost::shared_ptr<AerosolExtinction>& A, const std::string& Aname)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new AerosolParamOutput
     (boost::dynamic_pointer_cast<AerosolExtinctionImpBase>(A), Aname));
}
REGISTER_LUA_DERIVED_CLASS(AerosolParamOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &aco_create),
 luabind::def("create", &aco_create2)
]
REGISTER_LUA_END()
#endif

// See base class for description

void AerosolParamOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the state
  boost::shared_ptr<AerosolExtinctionImpBase> afreeze = 
    boost::dynamic_pointer_cast<AerosolExtinctionImpBase>(a->clone());
  std::string aer_name = aname;
  if(aer_name == "")
    aer_name = afreeze->aerosol_name();
  std::string mod_name = afreeze->model_short_name();
  boost::algorithm::to_lower(aer_name);
  out->register_data_source
    ("/RetrievalResults/aerosol_" + aer_name + "_" + mod_name + "_param_apriori",
     &AerosolExtinctionImpBase::aerosol_parameter, afreeze);
}

void AerosolParamOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  std::string aer_name = aname;
  if(aer_name == "")
    aer_name = a->aerosol_name();
  std::string mod_name = a->model_short_name();
  boost::algorithm::to_lower(aer_name);
  out->register_data_source
    ("/RetrievalResults/aerosol_" + aer_name + "_" + mod_name + "_param",
     &AerosolExtinctionImpBase::aerosol_parameter, a);
  out->register_data_source
    ("/RetrievalResults/aerosol_" + aer_name + "_" + mod_name + "_param_uncert",
    &AerosolExtinctionImpBase::aerosol_parameter_uncertainty, a);
}

