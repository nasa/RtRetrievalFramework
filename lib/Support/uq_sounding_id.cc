#include "uq_sounding_id.h"
#include "fp_exception.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(UqSoundingId, HdfSoundingId)
.def(luabind::constructor<const HdfFile&, const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Read Hdf file and determine sounding id information.
//-----------------------------------------------------------------------

UqSoundingId::UqSoundingId(const HdfFile& File, const std::string& Sounding_id)
{
    initialize(File, Sounding_id);
}

//-----------------------------------------------------------------------
/// Initialize sounding id information
//-----------------------------------------------------------------------

void UqSoundingId::initialize(const HdfFile& File, const std::string& Sounding_id)
{
    try {
        sounding_id_int = boost::lexical_cast<int64_t>(Sounding_id);

        Array<int64_t, 1> sid =
            File.read_field<int64_t, 1>("SoundingGeometry/sounding_id");

        //-----------------------------------------------------------------------
        // Search for sounding id.
        //-----------------------------------------------------------------------

        bool found = false;

        for(int i = 0; !found && i < sid.rows(); ++i) {
            if(sid(i) == sounding_id_int) {
                frame_number_ = i;
                found = true;
            }
        }

        // Sounding number is always 0, no dimension
        sounding_number_ = 0;

        if(!found) {
            Exception e;
            e << "The sounding " << Sounding_id << " was not found";
            throw e;
        }
    } catch(Exception& E) {
        E << " in the file " << File.file_name();
        throw E;
    }
}


void UqSoundingId::print(std::ostream& Os) const
{
    Os << "UqSoundingId: \n"
       << "  Frame number:    " << frame_number() << "\n"
       << "  Sounding ID:     " << sounding_id() << "\n";
}
