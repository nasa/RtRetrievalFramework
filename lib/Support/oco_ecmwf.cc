#include "oco_ecmwf.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(OcoEcmwf, Meteorology)
.def(luabind::constructor<std::string, const boost::shared_ptr<HdfSoundingId>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
/// \param Fname File to open
/// \param Hdf_sounding_id The sounding id to read in the file.
/// and pressure as the average value for all the sounding numbers.
//-----------------------------------------------------------------------

OcoEcmwf::OcoEcmwf(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& Hdf_sounding_id)
: h(Fname), hsid(Hdf_sounding_id)
{
}

//-----------------------------------------------------------------------
/// Add the capability to return Ozone VMR
//-----------------------------------------------------------------------

blitz::Array<double, 1> OcoEcmwf::vmr(const std::string& Species) const
{
    std::string species_upper = Species;
    boost::algorithm::to_upper(species_upper);
    if (species_upper == "O3") {
        return ozone_vmr();
    } else {
        return Meteorology::vmr(Species);
    }   
}

//-----------------------------------------------------------------------
/// Return the Ozone VMR.
//-----------------------------------------------------------------------

blitz::Array<double, 1> OcoEcmwf::ozone_vmr() const
{
    Array<double, 1> s = ozone_mmr();
    Array<double, 1> vmr(s.shape());
    vmr = s * OldConstant::molar_weight_dry_air / OldConstant::molar_weight_ozone;
    return vmr;
}

//-----------------------------------------------------------------------
/// Read a field where a single number is expected to be returned
//-----------------------------------------------------------------------

double OcoEcmwf::read_scalar(const std::string& Field) const
{
  return h.read_field<double, 2>
    (Field,
     TinyVector<int, 2>(hsid->frame_number(), hsid->sounding_number()),
     TinyVector<int, 2>(1,1))(0,0);
}

//-----------------------------------------------------------------------
/// Read a field and the pressure it is reported on. Average if needed.
//-----------------------------------------------------------------------

Array<double, 1> OcoEcmwf::read_array(const std::string& Field) const
{
  firstIndex i1; secondIndex i2;
  TinyVector<int, 3> sz = h.read_shape<3>(Field);
  Array<double, 3> traw = h.read_field<double, 3>
    (Field,
     TinyVector<int, 3>(hsid->frame_number(), hsid->sounding_number(), 0),
     TinyVector<int, 3>(1,1,sz[2]));
  Array<double, 1> V(traw.extent(thirdDim));

  V = traw(0, 0, Range::all());

  return V;
}
