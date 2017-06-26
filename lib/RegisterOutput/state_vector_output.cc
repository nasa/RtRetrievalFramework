#include "state_vector_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(StateVectorOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<StateVector>& >())
REGISTER_LUA_END()
#endif

// See base class for description
void StateVectorOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  blitz::Array<double, 1> initial_sv = sv->state().copy();
  out->register_data_source("/RetrievedStateVector/state_vector_initial", initial_sv);
}

void StateVectorOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source("/RetrievedStateVector/state_vector_names", 
	   &StateVector::state_vector_name, sv);
}

