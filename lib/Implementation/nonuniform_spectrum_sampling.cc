#include "nonuniform_spectrum_sampling.h"
#include "ostream_pad.h"
#include <algorithm>

using namespace FullPhysics;
using namespace std;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NonuniformSpectrumSampling, SpectrumSampling)
.def(luabind::constructor<const HeritageFile&,
                          const HeritageFile&,
                          const HeritageFile&,
                          const boost::shared_ptr<SpectrumSampling>&>())
.def(luabind::constructor<const std::string&,
                          const std::string&,
                          const std::string&,
                          const boost::shared_ptr<SpectrumSampling>&>())
.def(luabind::constructor<const SpectralDomain&,
                          const SpectralDomain&,
                          const SpectralDomain&,
                          const boost::shared_ptr<SpectrumSampling>&>())
REGISTER_LUA_END()
#endif



//-----------------------------------------------------------------------
/// Constructor. This creates a grid with no assumption made
/// on the uniformity of the spacing.
//-----------------------------------------------------------------------

NonuniformSpectrumSampling::NonuniformSpectrumSampling(
const SpectralDomain& Grid,
const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling
)
: SpectrumSampling(1), interpolated_sampling(Interpolated_sampling)
{
  spec_domain.push_back(sort_sd(Grid));
}

//-----------------------------------------------------------------------
/// Constructor for 3 spectrum (this is a useful case, because it
/// matches GOSAT)
//-----------------------------------------------------------------------

NonuniformSpectrumSampling::NonuniformSpectrumSampling(
const SpectralDomain& Grid1,
const SpectralDomain& Grid2,
const SpectralDomain& Grid3,
const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling
)
: SpectrumSampling(3), interpolated_sampling(Interpolated_sampling)
{
  spec_domain.push_back(sort_sd(Grid1));
  spec_domain.push_back(sort_sd(Grid2));
  spec_domain.push_back(sort_sd(Grid3));
}

//-----------------------------------------------------------------------
/// Constructor. This creates a grid with no assumption made
/// on the uniformity of the spacing. Grid points are read
/// from a text file with heritage-file format.
//-----------------------------------------------------------------------

NonuniformSpectrumSampling::NonuniformSpectrumSampling(
const HeritageFile& Grid_file,
const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling
)
: SpectrumSampling(1), interpolated_sampling(Interpolated_sampling)
{
  spec_domain.push_back(sort_sd(Grid_file.data("grid_points_wn")));
}

//-----------------------------------------------------------------------
/// Another constructor for 3 spectrum (this is a useful case, 
/// because it matches GOSAT)
//-----------------------------------------------------------------------

NonuniformSpectrumSampling::NonuniformSpectrumSampling(
const HeritageFile& Grid_file1,
const HeritageFile& Grid_file2,
const HeritageFile& Grid_file3,
const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling
)
: SpectrumSampling(3), interpolated_sampling(Interpolated_sampling)
{
  spec_domain.push_back(sort_sd(Grid_file1.data("grid_points_wn")));
  spec_domain.push_back(sort_sd(Grid_file2.data("grid_points_wn")));
  spec_domain.push_back(sort_sd(Grid_file3.data("grid_points_wn")));
}


//-----------------------------------------------------------------------
/// Constructor. This creates a grid with no assumption made
/// on the uniformity of the spacing. Grid points are read
/// from a text file with heritage-file format.
//-----------------------------------------------------------------------

NonuniformSpectrumSampling::NonuniformSpectrumSampling(
const std::string& Grid_file,
const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling
)
: SpectrumSampling(1), interpolated_sampling(Interpolated_sampling)
{
  if(Grid_file == "") {
    Array<double, 1> empty;
    spec_domain.push_back(SpectralDomain(empty));
  } else {
    HeritageFile f(Grid_file);
    spec_domain.push_back(sort_sd(f.data("grid_points_wn")));
  }
}

//-----------------------------------------------------------------------
/// Another constructor for 3 spectrum (this is a useful case, 
/// because it matches GOSAT)
//-----------------------------------------------------------------------

NonuniformSpectrumSampling::NonuniformSpectrumSampling(
const std::string& Grid_file1,
const std::string& Grid_file2,
const std::string& Grid_file3,
const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling
)
: SpectrumSampling(3), interpolated_sampling(Interpolated_sampling)
{
  if(Grid_file1 == "") {
    Array<double, 1> empty;
    spec_domain.push_back(SpectralDomain(empty));
  } else {
    HeritageFile f(Grid_file1);
    spec_domain.push_back(sort_sd(f.data("grid_points_wn")));
  }
  if(Grid_file2 == "") {
    Array<double, 1> empty;
    spec_domain.push_back(SpectralDomain(empty));
  } else {
    HeritageFile f(Grid_file2);
    spec_domain.push_back(sort_sd(f.data("grid_points_wn")));
  }
  if(Grid_file3 == "") {
    Array<double, 1> empty;
    spec_domain.push_back(SpectralDomain(empty));
  } else {
    HeritageFile f(Grid_file3);
    spec_domain.push_back(sort_sd(f.data("grid_points_wn")));
  }
}

SpectralDomain NonuniformSpectrumSampling::sort_sd
(const SpectralDomain& In) const
{
  Array<double, 1> res(In.data().copy());
  sort(res.data(), res.data()+res.size());
  return SpectralDomain(res, In.units());
}

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

void NonuniformSpectrumSampling::print(std::ostream& Os) const 
{ 
  OstreamPad opad(Os, "    ");
  Os << "NonuniformSpectrumSampling\n";
  Os << "  Uniform spectrum sampling:\n";
  opad << *interpolated_sampling << "\n";
  opad.strict_sync();
  for(int i = 0; i < number_spectrometer(); ++i)
    if(spec_domain[i].data().rows() <= 0)
      Os << "  Band " << i + 1 << ":\n"
         << "     Nonuniform RT grid not provided, using uniform grid for RT\n";
    else
      Os << "  Band " << i + 1 << ":\n"
         << "     grid_start: " << spec_domain[i].data()(0) << "\n"
         << "     grid_end:   " << spec_domain[i].data()(spec_domain[i].data().rows() - 1) << "\n"
         << "     grid_points: " << spec_domain[i].data().rows() << "\n";
}

// See base class for description
SpectralDomain NonuniformSpectrumSampling::spectral_domain
(int spec_index,
 const SpectralDomain& Lowres_grid, 
 const DoubleWithUnit& Ils_half_width) const
{
  range_check(spec_index, 0, number_spectrometer());
  SpectralDomain ispec = interpolated_sampling->
    spectral_domain(spec_index, Lowres_grid, Ils_half_width);
  // If we aren't actually doing nonuniform sampling, just return the
  // spectral domain we want to sample to.
  if(spec_domain[spec_index].data().rows() == 0)
    return ispec;
  // Find all points in ispec that are closest to one of the grid points.
  std::map<double, bool> res;
  typedef std::map<double, bool>::iterator it_type;
  typedef std::map<double, bool>::value_type val_type;
  BOOST_FOREACH(double x, ispec.data())
    res[x] = false;
  BOOST_FOREACH(double x, spec_domain[spec_index].convert_wave(ispec.units())) {
    it_type i = res.lower_bound(x);
    if(i != res.begin()) {
      it_type i2 = i;
      --i2;
      if(fabs(x - (*i2).first) < fabs(x - (*i).first))
	i = i2;
    }
    (*i).second = true;
  }
  // And use those points to make up the spectral domain.
  std::vector<double> resx;
  BOOST_FOREACH(const val_type& i, res)
    if(i.second)
      resx.push_back(i.first);
  Array<double, 1> sd(&resx[0], shape((int) resx.size()), duplicateData);
  return SpectralDomain(sd, ispec.units());
}
