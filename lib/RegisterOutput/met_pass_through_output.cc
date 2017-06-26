#include "met_pass_through_output.h"
#include "fp_exception.h"

// This won't be needed once version 3 becomes the default for boost filesystem
#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(MetPassThroughOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&>())
REGISTER_LUA_END()
#endif

void MetPassThroughOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    out->register_data_source
        ("/RetrievalResults/wind_speed_u_met", &Meteorology::windspeed_u, met_);

    out->register_data_source
        ("/RetrievalResults/wind_speed_v_met", &Meteorology::windspeed_v, met_);
}
