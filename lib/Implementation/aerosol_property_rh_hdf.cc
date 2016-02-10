#include "aerosol_property_rh_hdf.h"
#include "ostream_pad.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolPropertyRhHdf, AerosolProperty)
.def(luabind::constructor<const HdfFile&, const std::string&, 
     const boost::shared_ptr<Pressure>&, 
     const boost::shared_ptr<RelativeHumidity>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Read the given group in the given file for the aerosol properties.
//-----------------------------------------------------------------------

AerosolPropertyRhHdf::AerosolPropertyRhHdf
(const HdfFile& F, 
 const std::string& Group_name,
 const boost::shared_ptr<Pressure>& Press,
 const boost::shared_ptr<RelativeHumidity>& Rh
)
: rh(Rh), hdf_file(F.file_name()), hdf_group(Group_name)
{
  press = Press;
  Array<double, 1> wn(F.read_field<double, 1>(Group_name + "/wave_number"));
  Array<double, 1> rh(F.read_field<double, 1>(Group_name + "/relative_humidity"));
  Array<double, 2> 
    qscatv(F.read_field<double, 2>(Group_name + "/scattering_coefficient"));
  Array<double, 2> 
    qextv(F.read_field<double, 2>(Group_name + "/extinction_coefficient"));
  Array<double, 4>
    pfv(F.read_field<double, 4>(Group_name + "/phase_function_moment"));
  Array<double, 1> qscatvt(qscatv.rows());
  Array<double, 1> qextvt(qextv.rows());
  qscatvt = qscatv(Range::all(), 0);
  qextvt = qextv(Range::all(), 0);
  qext.reset(new LinearInterpolate<double, double>(wn.begin(), wn.end(), 
					   qextvt.begin()));
  qscat.reset(new LinearInterpolate<double, double>(wn.begin(), wn.end(), 
					    qscatvt.begin()));
  std::vector<Array<double, 2> > pf_vec;
  for(int i = 0; i < pfv.rows(); ++i)
    pf_vec.push_back(Array<double, 2>(pfv(i, 0, Range::all(), Range::all())));
  pf.reset(new ScatteringMomentInterpolate(wn.begin(), wn.end(),
					   pf_vec.begin()));
}

boost::shared_ptr<AerosolProperty> AerosolPropertyRhHdf::clone() const
{
  return clone(press->clone(), rh->clone());
}

boost::shared_ptr<AerosolProperty> AerosolPropertyRhHdf::clone
(const boost::shared_ptr<Pressure>& Press,
 const boost::shared_ptr<RelativeHumidity>& Rh) const
{
  HdfFile f(hdf_file);
  return boost::shared_ptr<AerosolProperty>
    (new AerosolPropertyRhHdf(f, hdf_group, Press, Rh));
}

ArrayAd<double, 1> AerosolPropertyRhHdf::extinction_coefficient_each_layer
(double wn) const
{
  firstIndex i1; secondIndex i2;
  AutoDerivative<double> t = (*qext)(wn);
  ArrayAd<double, 1> res(press->number_layer(), t.number_variable());
  res.value() = t.value();
  if(t.number_variable() > 0)
    res.jacobian() = t.gradient()(i2);
  return res;
}

ArrayAd<double, 1> AerosolPropertyRhHdf::scattering_coefficient_each_layer
(double wn) const
{
  firstIndex i1; secondIndex i2;
  AutoDerivative<double> t = (*qscat)(wn);
  ArrayAd<double, 1> res(press->number_layer(), t.number_variable());
  res.value() = t.value();
  if(t.number_variable() > 0)
    res.jacobian() = t.gradient()(i2);
  return res;
}

ArrayAd<double, 3> AerosolPropertyRhHdf::phase_function_moment_each_layer
(double wn, int nmom, int nscatt) const
{ 
  firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
  ArrayAd<double, 2> t = (*pf)(wn, nmom, nscatt);
  ArrayAd<double, 3> res(t.rows(), press->number_layer(), 
			 t.cols(), t.number_variable());
  res.value() = t.value()(i1, i3);
  if(t.number_variable() > 0)
    res.jacobian() = t.jacobian()(i1, i3, i4);
  return res; 
}

void AerosolPropertyRhHdf::print(std::ostream& Os) const 
{ 
  Os << "AerosolPropertyRhHdf:\n"
     << "  Hdf file:  " << hdf_file << "\n"
     << "  Hdf group: " << hdf_group << "\n";
}
