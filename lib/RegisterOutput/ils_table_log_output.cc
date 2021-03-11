#include "ils_table_log_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> ils_table_create_log
(const boost::shared_ptr<IlsFunction>& Ils, const std::string& Hdf_band_name)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new IlsTableLogOutput
     (boost::dynamic_pointer_cast<IlsTableLog>(Ils), Hdf_band_name));
}
REGISTER_LUA_DERIVED_CLASS(IlsTableLogOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<IlsTableLog>&,
			  std::string>())
.scope
[
 luabind::def("create", &ils_table_create_log)
]
REGISTER_LUA_END()
#endif

// See base class for description

void IlsTableLogOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the instrument state
  boost::shared_ptr<IlsTableLog> dfreeze =
    boost::dynamic_pointer_cast<IlsTableLog>(d->clone());
  out->register_data_source
    ("/RetrievalResults/ils_scale_apriori_" + hdf_band_name,
     &IlsTableLog::ils_scale, dfreeze);
}

void IlsTableLogOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source
    ("/RetrievalResults/ils_scale_" + hdf_band_name,
     &IlsTableLog::ils_scale, d);
  out->register_data_source
    ("/RetrievalResults/ils_scale_uncert_" + hdf_band_name,
     &IlsTableLog::ils_scale_uncertainty, d);
}
