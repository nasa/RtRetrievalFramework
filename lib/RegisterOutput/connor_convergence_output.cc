#include "connor_convergence_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> conn_conv_create
(const boost::shared_ptr<ConvergenceCheck>& C)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new ConnorConvergenceOutput
     (boost::dynamic_pointer_cast<ConnorConvergence>(C)));
}
REGISTER_LUA_DERIVED_CLASS(ConnorConvergenceOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<ConnorConvergence>& >())
.scope
[
 luabind::def("create", &conn_conv_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void ConnorConvergenceOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source("/Metadata/RetrievalIterationLimit", 
			   &ConnorConvergence::maximum_number_iteration, 
			   conv);
}

