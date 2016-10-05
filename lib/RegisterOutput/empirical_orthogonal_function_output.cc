#include "empirical_orthogonal_function_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> eof_create
(const boost::shared_ptr<InstrumentCorrection>& Ic, 
 const std::string& Hdf_band_name)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new EmpiricalOrthogonalFunctionOutput
     (boost::dynamic_pointer_cast<EmpiricalOrthogonalFunction>(Ic), 
      Hdf_band_name));
}
REGISTER_LUA_DERIVED_CLASS(EmpiricalOrthogonalFunctionOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &eof_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void EmpiricalOrthogonalFunctionOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the instrument state
  boost::shared_ptr<EmpiricalOrthogonalFunction> eoffreeze = 
    boost::dynamic_pointer_cast<EmpiricalOrthogonalFunction>(eof->clone());
  out->register_data_source
    ("/RetrievalResults/eof_" +
     boost::lexical_cast<std::string>(eof->order()) +
     "_scale_apriori_" + hdf_band_name, 
     &EmpiricalOrthogonalFunction::scale, eoffreeze);
}

void EmpiricalOrthogonalFunctionOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source
    ("/RetrievalResults/eof_" +
     boost::lexical_cast<std::string>(eof->order()) +
     "_scale_" + hdf_band_name, 
     &EmpiricalOrthogonalFunction::scale, eof);
  out->register_data_source
    ("/RetrievalResults/eof_" +
     boost::lexical_cast<std::string>(eof->order()) +
     "_scale_uncert_" + hdf_band_name, 
     &EmpiricalOrthogonalFunction::scale_uncertainty, eof);
}

