#include "acos_met_file.h"
#include "acos_sounding_id.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AcosMetFile, Meteorology)
.def(luabind::constructor<std::string, 
			  const boost::shared_ptr<HdfSoundingId>&,
			  bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
/// \param Fname File to open
/// \param Hdf_sounding_id The sounding id to read in the file.
/// \param Avg_sounding_number If true, then we get the temperature
/// and pressure as the average value for all the sounding numbers.
//-----------------------------------------------------------------------

AcosMetFile::AcosMetFile(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& 
		     Hdf_sounding_id, bool Avg_sounding_number)
: h(Fname), hsid(Hdf_sounding_id), average_sounding_number(Avg_sounding_number)
{
  if (h.has_object("/Meteorology/specific_humidity_profile_met"))
    met_format = true;
  else
    met_format = false;
}

//-----------------------------------------------------------------------
/// \deprecated
/// This constructor is used by the deprecated class
/// InitialGuessHeritage only. Once that class has been removed we can
/// get rid of this constructor.
///
/// Constructor that reads the sounding information from the run
/// config file
//-----------------------------------------------------------------------

AcosMetFile::AcosMetFile(const std::string& Fname, const HeritageFile& Run_file)
  : h(Fname), met_format(false)
{
  std::string sid = Run_file.value<std::string>("SOUNDING_INFO/sounding_id");
  HdfFile h(Run_file.file_value("SOUNDING_INFO/spectrum_file"));
  std::vector<boost::shared_ptr<HdfSoundingId> > sidv = 
    AcosSoundingId::create(h, sid);
  average_sounding_number = (sidv.size() > 1);
  hsid = sidv[0];
}

//-----------------------------------------------------------------------
/// Read a field where a single number is expected to be returned
//-----------------------------------------------------------------------

double AcosMetFile::read_scalar(const std::string& Field) const
{
  int spec_index = 0;
  TinyVector<int, 3> sz = h.read_shape<3>(Field);
  Array<double, 3> data = h.read_field<double, 3>
    (Field,
     TinyVector<int, 3>(hsid->frame_number(), spec_index, 0),
     TinyVector<int, 3>(1,1,sz[2]));
  if(average_sounding_number)
    return sum(data(0, 0, Range::all())) / 2;
  else
    return data(0, 0, hsid->sounding_number());
}

//-----------------------------------------------------------------------
/// Read a field and the pressure it is reported on. Average if needed.
//-----------------------------------------------------------------------

blitz::Array<double, 1> AcosMetFile::read_array(const std::string& Field) const
{
  firstIndex i1; secondIndex i2;
  int spec_index = 0;
  TinyVector<int, 4> sz = h.read_shape<4>(Field);
  Array<double, 4> traw = h.read_field<double, 4>
    (Field,
     TinyVector<int, 4>(hsid->frame_number(), spec_index, 0, 0),
     TinyVector<int, 4>(1,1,sz[2],sz[3]));
  Array<double, 1> V(traw.extent(fourthDim));
  if(average_sounding_number) {
    V = sum(traw(0, 0, Range::all(), Range::all())(i2, i1), i2) / 2;
  } else {
    V = traw(0, 0, hsid->sounding_number(), Range::all());
  }
  return V;
}
