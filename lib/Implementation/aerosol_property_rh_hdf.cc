#include "aerosol_property_rh_hdf.h"
#include "ostream_pad.h"
#include <boost/make_shared.hpp>
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
  Array<double, 1> rhv(F.read_field<double, 1>(Group_name + "/relative_humidity"));
  Array<double, 2> 
    qscatv(F.read_field<double, 2>(Group_name + "/scattering_coefficient"));
  Array<double, 2> 
    qextv(F.read_field<double, 2>(Group_name + "/extinction_coefficient"));
  Array<double, 4>
    pfv(F.read_field<double, 4>(Group_name + "/phase_function_moment"));
  if(rhv.rows() != qscatv.cols() ||
     rhv.rows() != qextv.cols() ||
     rhv.rows() != pfv.cols())
    throw Exception("Mismatch between the number of relative humidity values and the aerosol property arrays");
  Array<double, 1> qscatvt(qscatv.rows());
  Array<double, 1> qextvt(qextv.rows());
  for(int i = 0; i < rhv.rows(); ++i) {
    rh_val.push_back(rhv(i));
    rh_val_d.push_back(rhv(i));
    qscatvt = qscatv(Range::all(), i);
    qextvt = qextv(Range::all(), i);
    qext.push_back(boost::make_shared<LinearInterpolate<double, double> >
		  (wn.begin(), wn.end(), qextvt.begin()));
    qscat.push_back(boost::make_shared<LinearInterpolate<double, double> >
		   (wn.begin(), wn.end(), qscatvt.begin()));
    std::vector<Array<double, 2> > pf_vec;
    for(int j = 0; j < pfv.rows(); ++j)
      pf_vec.push_back(Array<double, 2>(pfv(j, i, Range::all(), Range::all())));
    pf.push_back(boost::make_shared<ScatteringMomentInterpolate>
		(wn.begin(), wn.end(), pf_vec.begin()));
  }
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

// Right now we don't support extinction_coefficient_each_layer 
// having a jacobian,
// this has to do with our optimization in the LIDORT jacobian
// calculation that uses the "intermediate" variables - see
// rt_atmosphere.h for more details on this. So we will just leave
// this jacobian out. Not clear what if any impact this will have on
// the retrieval, but we'll do this to try things out.

ArrayAd<double, 1> AerosolPropertyRhHdf::extinction_coefficient_each_layer_not_used
(double wn) const
{
  std::vector<AutoDerivative<double> > qextv;
  for(int i = 0; i < (int) rh_val.size(); ++i)
    qextv.push_back((*qext[i])(wn));
  LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > 
    lin(rh_val.begin(), rh_val.end(), qextv.begin());
  ArrayAd<double, 1> rhl = rh->relative_humidity_layer();
  blitz::Array<AutoDerivative<double>, 1> res(rhl.rows());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = lin(rhl(i));
  return ArrayAd<double, 1>(res);
}

ArrayAd<double, 1> AerosolPropertyRhHdf::extinction_coefficient_each_layer
(double wn) const
{
  std::vector<double> qextv;
  for(int i = 0; i < (int) rh_val.size(); ++i)
    qextv.push_back((*qext[i])(wn));
  LinearInterpolate<double, double> 
    lin(rh_val_d.begin(), rh_val_d.end(), qextv.begin());
  blitz::Array<double, 1> rhl = rh->relative_humidity_layer().value();
  blitz::Array<double, 1> res(rhl.rows());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = lin(rhl(i));
  return ArrayAd<double, 1>(res);
}

// Right now we don't support scattering_coefficient_each_layer 
// having a jacobian,
// this has to do with our optimization in the LIDORT jacobian
// calculation that uses the "intermediate" variables - see
// rt_atmosphere.h for more details on this. So we will just leave
// this jacobian out. Not clear what if any impact this will have on
// the retrieval, but we'll do this to try things out.

ArrayAd<double, 1> AerosolPropertyRhHdf::scattering_coefficient_each_layer_not_used
(double wn) const
{
  std::vector<AutoDerivative<double> > qscatv;
  for(int i = 0; i < (int) rh_val.size(); ++i)
    qscatv.push_back((*qscat[i])(wn));
  LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > 
    lin(rh_val.begin(), rh_val.end(), qscatv.begin());
  ArrayAd<double, 1> rhl = rh->relative_humidity_layer();
  blitz::Array<AutoDerivative<double>, 1> res(rhl.rows());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = lin(rhl(i));
  return ArrayAd<double, 1>(res);
}

ArrayAd<double, 1> AerosolPropertyRhHdf::scattering_coefficient_each_layer
(double wn) const
{
  std::vector<double> qscatv;
  for(int i = 0; i < (int) rh_val.size(); ++i)
    qscatv.push_back((*qscat[i])(wn));
  LinearInterpolate<double, double> 
    lin(rh_val_d.begin(), rh_val_d.end(), qscatv.begin());
  blitz::Array<double, 1> rhl = rh->relative_humidity_layer().value();
  blitz::Array<double, 1> res(rhl.rows());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = lin(rhl(i));
  return ArrayAd<double, 1>(res);
}

// Right now we don't support phase_function_moment having a jacobian,
// this has to do with our optimization in the LIDORT jacobian
// calculation that uses the "intermediate" variables - see
// rt_atmosphere.h for more details on this. So we will just leave
// this jacobian out. Not clear what if any impact this will have on
// the retrieval, but we'll do this to try things out.
ArrayAd<double, 3> AerosolPropertyRhHdf::phase_function_moment_each_layer_not_used
(double wn, int nmom, int nscatt) const
{ 
  firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
  std::vector<blitz::Array<AutoDerivative<double>, 2> > pfv;
  for(int i = 0; i < (int) rh_val.size(); ++i) {
    blitz::Array<double, 2> t((*pf[i])(wn, nmom, nscatt));
    blitz::Array<AutoDerivative<double>, 2> t2(t.shape());
    for(int j  = 0; j < t2.rows(); ++j)
      for(int k  = 0; k < t2.cols(); ++k)
	t2(j, k) = t(j, k);
    pfv.push_back(t2);
  }
  LinearInterpolate<AutoDerivative<double>, 
    blitz::Array<AutoDerivative<double>, 2> >
    lin(rh_val.begin(), rh_val.end(), pfv.begin());
  ArrayAd<double, 1> rhl = rh->relative_humidity_layer();
  blitz::Array<AutoDerivative<double>, 3> 
    res(pfv[0].rows(), press->number_layer(), pfv[1].cols());
  for(int i = 0; i < res.cols(); ++i)
    res(Range::all(), i, Range::all()) = lin(rhl(i));
  return ArrayAd<double, 3>(res); 
}

ArrayAd<double, 3> AerosolPropertyRhHdf::phase_function_moment_each_layer
(double wn, int nmom, int nscatt) const
{ 
  firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
  std::vector<blitz::Array<double, 2> > pfv;
  for(int i = 0; i < (int) rh_val.size(); ++i) {
    blitz::Array<double, 2> t((*pf[i])(wn, nmom, nscatt));
    pfv.push_back(t);
  }
  LinearInterpolate<double, blitz::Array<double, 2> >
    lin(rh_val_d.begin(), rh_val_d.end(), pfv.begin());
  blitz::Array<double, 1> rhl = rh->relative_humidity_layer().value();
  blitz::Array<double, 3> 
    res(pfv[0].rows(), press->number_layer(), pfv[1].cols());
  for(int i = 0; i < res.cols(); ++i)
    res(Range::all(), i, Range::all()) = lin(rhl(i));
  return ArrayAd<double, 3>(res); 
}

void AerosolPropertyRhHdf::print(std::ostream& Os) const 
{ 
  Os << "AerosolPropertyRhHdf:\n"
     << "  Hdf file:  " << hdf_file << "\n"
     << "  Hdf group: " << hdf_group << "\n";
}
