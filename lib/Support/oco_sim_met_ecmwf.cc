#include "oco_sim_met_ecmwf.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(OcoSimMetEcmwf, Ecmwf)
.def(luabind::constructor<std::string, 
			  const boost::shared_ptr<HdfSoundingId>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
/// \param Fname File to open
/// \param Hdf_sounding_id The sounding id to read in the file.
/// and pressure as the average value for all the sounding numbers.
//-----------------------------------------------------------------------

OcoSimMetEcmwf::OcoSimMetEcmwf(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& Hdf_sounding_id)
: h(Fname), hsid(Hdf_sounding_id)
{
}

//-----------------------------------------------------------------------
/// Read a field where a single number is expected to be returned
//-----------------------------------------------------------------------

double OcoSimMetEcmwf::read(const std::string& Field) const
{
  return h.read_field<double, 2>
    (Field,
     TinyVector<int, 2>(hsid->frame_number(), hsid->sounding_number()),
     TinyVector<int, 2>(1,1))(0,0);
}

//-----------------------------------------------------------------------
/// Read a field and the pressure it is reported on. Average if needed.
//-----------------------------------------------------------------------

void OcoSimMetEcmwf::read(const std::string& Field, blitz::Array<double, 1>& P, 
			  blitz::Array<double, 1>& V) const
{
  firstIndex i1; secondIndex i2;
  TinyVector<int, 3> sz = h.read_shape<3>(Field);
  Array<double, 3> traw = h.read_field<double, 3>
    (Field,
     TinyVector<int, 3>(hsid->frame_number(), hsid->sounding_number(), 0),
     TinyVector<int, 3>(1,1,sz[2]));
  TinyVector<int, 3> sz2 = h.read_shape<3>(Field + "_pressures");
  Array<double, 3> praw = h.read_field<double, 3>
    (Field + "_pressures",
     TinyVector<int, 3>(hsid->frame_number(), hsid->sounding_number(), 0),
     TinyVector<int, 3>(1,1,sz2[2]));

  V.resize(traw.extent(thirdDim));
  P.resize(traw.extent(thirdDim));

  V = traw(0, 0, Range::all());
  P = praw(0, 0, Range::all());
}
