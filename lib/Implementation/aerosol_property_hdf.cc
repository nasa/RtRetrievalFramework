#include "aerosol_property_hdf.h"
#include "ostream_pad.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolPropertyHdf, AerosolProperty)
.def(luabind::constructor<const HdfFile&, const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Read the given group in the given file for the aerosol properties.
//-----------------------------------------------------------------------

AerosolPropertyHdf::AerosolPropertyHdf(const HdfFile& F, 
				       const std::string& Group_name)
: hdf_file(F.file_name()), hdf_group(Group_name)
{
  Array<double, 1> wn(F.read_field<double, 1>(Group_name + "/wave_number"));
  Array<double, 1> 
    qscatv(F.read_field<double, 1>(Group_name + "/scattering_coefficient"));
  Array<double, 1> 
    qextv(F.read_field<double, 1>(Group_name + "/extinction_coefficient"));
  Array<double, 3>
    pfv(F.read_field<double, 3>(Group_name + "/phase_function_moment"));
  qext.reset(new LinearInterpolate<double, double>(wn.begin(), wn.end(), 
					   qextv.begin()));
  qscat.reset(new LinearInterpolate<double, double>(wn.begin(), wn.end(), 
					    qscatv.begin()));
  std::vector<Array<double, 2> > pf_vec;
  for(int i = 0; i < pfv.rows(); ++i)
    pf_vec.push_back(Array<double, 2>(pfv(i, Range::all(), Range::all())));
  pf.reset(new ScatteringMomentInterpolate(wn.begin(), wn.end(),
					   pf_vec.begin()));
}

void AerosolPropertyHdf::print(std::ostream& Os) const 
{ 
  Os << "AerosolPropertyHdf:\n"
     << "  Hdf file:  " << hdf_file << "\n"
     << "  Hdf group: " << hdf_group << "\n";
}
