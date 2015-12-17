#include "oco_sounding_id.h"
#include "fp_exception.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(OcoSoundingId,HdfSoundingId)
.def(luabind::constructor<const HdfFile&, const std::string&>())
REGISTER_LUA_END()
#endif


//-----------------------------------------------------------------------
/// Read Hdf file and determine sounding id information.
///
/// \param File The file to read.
/// \param Sounding_id  The sounding ID.
///
/// The sounding ID should be in the format:
///  -# 200907261717015 - Where the ending 5 is the sounding number
//-----------------------------------------------------------------------

OcoSoundingId::OcoSoundingId(const HdfFile& File, 
			     const std::string& Sounding_id)
{
  initialize(File, Sounding_id);
}

//-----------------------------------------------------------------------
/// Do the actual work.
//-----------------------------------------------------------------------

void OcoSoundingId::initialize(const HdfFile& File, 
			       const std::string& Sounding_id)
{
  try {
    sounding_id_int = boost::lexical_cast<int64_t>(Sounding_id);

    Array<int64_t, 2> sid = 
      File.read_field<int64_t, 2>("SoundingGeometry/sounding_id");

    number_sounding = sid.cols();
    
    //-----------------------------------------------------------------------
    // Search for sounding id.
    //-----------------------------------------------------------------------

    bool found = false;
    for(int i = 0; !found && i < sid.rows(); ++i)
      for(int j = 0; !found && j < sid.cols(); ++j)
	if(sid(i, j) == sounding_id_int) {
	  frame_number_ = i;
	  sounding_number_ = j;
	  found = true;
	}
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


void OcoSoundingId::print(std::ostream& Os) const
{
  Os << "OcoSoundingId: \n" 
     << "  Frame number:    " << frame_number() << "\n"
     << "  Sounding number: " << sounding_number() << "\n"
     << "  Sounding ID:     " << sounding_id() << "\n";
}
