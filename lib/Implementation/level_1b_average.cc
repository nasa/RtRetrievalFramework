#include "level_1b_average.h"
#include <boost/foreach.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(Level1bAverage, Level1b)
.def(luabind::constructor<const std::vector<boost::shared_ptr<Level1b> > &>())
REGISTER_LUA_END()
#endif

DoubleWithUnit Level1bAverage::latitude(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->latitude(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->latitude(i);
  return sum / ((int) l1b.size());
}

DoubleWithUnit Level1bAverage::longitude(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->longitude(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->longitude(i);
  return sum / ((int) l1b.size());
}

DoubleWithUnit Level1bAverage::sounding_zenith(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->sounding_zenith(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->sounding_zenith(i);
  return sum / ((int) l1b.size());
}

DoubleWithUnit Level1bAverage::sounding_azimuth(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->sounding_azimuth(1).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->sounding_azimuth(i);
  return sum / ((int) l1b.size());
}

Array<double, 1> Level1bAverage::stokes_coefficient(int i) const
{
  Array<double, 1> res(l1b[0]->stokes_coefficient(i).rows());
  res = 0;
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    res += f->stokes_coefficient(i);
  res /= (int) l1b.size();
  return res;
}

DoubleWithUnit Level1bAverage::solar_zenith(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->solar_zenith(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->solar_zenith(i);
  return sum / ((int) l1b.size());
}

DoubleWithUnit Level1bAverage::solar_azimuth(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->solar_azimuth(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->solar_azimuth(i);
  return sum / ((int) l1b.size());
}

DoubleWithUnit Level1bAverage::altitude(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->altitude(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->altitude(i);
  return sum / ((int) l1b.size());
}

ArrayWithUnit<double, 1> Level1bAverage::spectral_coefficient(int i) const
{
  Array<double, 1> sum(l1b[0]->spectral_coefficient(i).value.rows());
  sum = 0.0;
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->spectral_coefficient(i).value;
  ArrayWithUnit<double, 1> res;
  res.value.resize(sum.rows());
  res.value = sum / ((int) l1b.size());
  res.units = l1b[0]->spectral_coefficient(i).units;
  return res;
}

SpectralRange Level1bAverage::radiance(int Spec_index) const
{
  SpectralRange t = l1b[0]->radiance(Spec_index);
  Array<double, 1> rad(t.data());
  Array<double, 1> sum(rad.shape());
  sum = 0;
  bool have_uncertainty = (t.uncertainty().rows() > 0);
  if(have_uncertainty)
    sum += sqr(t.uncertainty());
  for(int i = 1; i < (int) l1b.size(); ++i) {
    SpectralRange t2 = l1b[i]->radiance(Spec_index);
    rad += t2.data() * FullPhysics::conversion(t2.units(), t.units());
    if(have_uncertainty)
      sum += sqr(t2.uncertainty() * FullPhysics::conversion(t2.units(), t.units()));
  }
  rad /= ((int) l1b.size());
  Array<double, 1> uncer;
  if(have_uncertainty)
    uncer.reference(Array<double,1>(sqrt(sum) / (int) l1b.size()));
  return SpectralRange(rad, t.units(), uncer);
}

void Level1bAverage::print(std::ostream& Os) const
{
  Os << "Level1bAverage";
}
