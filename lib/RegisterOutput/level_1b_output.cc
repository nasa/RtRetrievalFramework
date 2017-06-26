#include "level_1b_output.h"

using namespace FullPhysics;


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(Level1bOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<Level1b>& >())
REGISTER_LUA_END()
#endif

// See base class for description

void Level1bOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source("/RetrievalHeader/exposure_index", 
			   &Level1b::exposure_index, f);
  out->register_data_source("/RetrievalHeader/sounding_id_reference", 
  			   &Level1b::sounding_id, f);
}

