#include "acos_sounding_id.h"
#include "fp_exception.h"
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/case_conv.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AcosSoundingId,HdfSoundingId)
.enum_("SoundingType")
[
 luabind::value("S_SOUNDING", AcosSoundingId::S_SOUNDING),
 luabind::value("P_SOUNDING", AcosSoundingId::P_SOUNDING)
]
.def(luabind::constructor<const HdfFile&, const std::string&, 
			  AcosSoundingId::SoundingType>())
.scope
[
 luabind::def("create", &AcosSoundingId::create)
]
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Read Hdf file and determine sounding id information.
///
/// \param File The file to read.
/// \param Sounding_id  The sounding ID.
/// \param Sounding_type Specify if we have S, P.
//-----------------------------------------------------------------------

AcosSoundingId::AcosSoundingId(const HdfFile& File, 
			       const std::string& Sounding_id, 
			       SoundingType Sounding_type)
{
  initialize(File, Sounding_id, Sounding_type);
}

//-----------------------------------------------------------------------
/// Parse a sounding ID string. This string indicates if we are doing
/// just a P retrieval, just a S retrieval, or averaging the two.
///
/// The sounding ID should be in the format:
///  -# 20090726171701S - use the S band
///  -# 20090726171701P - use the P band
///  -# 20090726171701  - Average the two bands
///
/// This return a vector of AcosSoundingId, length 1 if we are not averaging
/// or length 2 if we are averaging.
//-----------------------------------------------------------------------

std::vector<boost::shared_ptr<HdfSoundingId> >
AcosSoundingId::create(const HdfFile& File, 
			const std::string& Sounding_id)
{
  std::vector<boost::shared_ptr<HdfSoundingId> > res;
  boost::smatch pol_match;
  if(!boost::regex_match(Sounding_id, pol_match, 
			 boost::regex(" *(\\d+)([pPsS]?) *")))
    throw Exception("Unrecognized format sounding id " + Sounding_id);
  std::string sounding_id_base = pol_match[1];
  std::string pol_code = pol_match[2];
  boost::to_lower(pol_code);
  if(pol_code == "p")
    res.push_back(boost::shared_ptr<HdfSoundingId>
		  (new AcosSoundingId(File, sounding_id_base, P_SOUNDING)));
  else if(pol_code == "s")
    res.push_back(boost::shared_ptr<HdfSoundingId>
		  (new AcosSoundingId(File, sounding_id_base, S_SOUNDING)));
  else {
    res.push_back(boost::shared_ptr<HdfSoundingId>
		  (new AcosSoundingId(File, sounding_id_base, P_SOUNDING)));
    res.push_back(boost::shared_ptr<HdfSoundingId>
		  (new AcosSoundingId(File, sounding_id_base, S_SOUNDING)));
  }
  return res;
}

//-----------------------------------------------------------------------
/// Do the actual work.
//-----------------------------------------------------------------------

void AcosSoundingId::initialize(const HdfFile& File, 
			       const std::string& Sounding_id,
			       SoundingType Sounding_type)
{
  try {
    sounding_id_int = boost::lexical_cast<int64_t>(Sounding_id);
    switch(Sounding_type) {
    case S_SOUNDING:
      sounding_number_ = 1;
      break;
    case P_SOUNDING:
      sounding_id_int = boost::lexical_cast<int64_t>(Sounding_id);
      sounding_number_ = 0;
      break;
    default:
      throw Exception("Unrecognized Sounding_type");
    }

//-----------------------------------------------------------------------
/// Now, look up the sounding ID in the file to determine the frame number
//-----------------------------------------------------------------------

    Array<int64_t, 1> sid = 
      File.read_field<int64_t, 1>("SoundingHeader/sounding_id");
    for(frame_number_ = 0; sid(frame_number_) != sounding_id_int && 
	  frame_number_ < sid.extent(firstDim); ++frame_number_)
      ;
    if(frame_number_ ==sid.extent(firstDim)) {
      Exception e;
      e << "The sounding " << Sounding_id << " was not found";
      throw e;
    }
  } catch(Exception& E) {
    E << " in the file " << File.file_name();
    throw E;
  }
}


void AcosSoundingId::print(std::ostream& Os) const
{
  Os << "AcosSoundingId: \n" 
     << "  Frame number:    " << frame_number() << "\n"
     << "  Sounding number: " << sounding_number() << "\n"
     << "  Sounding ID:     " << sounding_id() << "\n";
}
