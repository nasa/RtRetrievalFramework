#include "spectral_window_range.h"
#include "dispersion.h"
#include "fp_exception.h"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
typedef void (SpectralWindowRange::*cfunc)(const std::vector<boost::shared_ptr<Dispersion> >&);
typedef const ArrayWithUnit<double, 3>& (SpectralWindowRange::*a1)(void) const;
REGISTER_LUA_DERIVED_CLASS(SpectralWindowRange, SpectralWindow)
.def(luabind::constructor<const ArrayWithUnit<double, 3>&>())
.def(luabind::constructor<const ArrayWithUnit<double, 3>&, const Array<bool, 2>&>())
.def(luabind::constructor<const ArrayWithUnit<double, 3>&, const Array<double, 2>&>())
.def("range_array", ((a1) &SpectralWindowRange::range_array))
.def("dispersion", ((cfunc) &SpectralWindowRange::dispersion))
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Construct a new window from the supplied information. We take
/// an array, which is number_spectrometer x number_microwindow x
/// 2. For a particular microwindow, we have 2 values, a lower and
/// upper range.
///
/// Note that we require the number of microwindows to be same for all
/// the spectrometers. This is just for convenience, it makes a
/// simpler interface. It is perfectly ok to have microwindow
/// with ranges of 0 to 0 - so you can simple set the number of
/// microwindows to whatever the maximum number is and use null range
/// microwindows for microwindows that aren't needed.
//-----------------------------------------------------------------------

SpectralWindowRange::SpectralWindowRange
(const ArrayWithUnit<double, 3>& Microwindow_ranges)
    : range_(Microwindow_ranges)
{
  if(range_.value.depth() != 2) {
    Exception e;
    e << "Microwindow_ranges needs to have a depth of 2\n"
      << "(for lower and upper bound). Found depth of " << range_.value.depth();
    throw e;
  }
  
}

//-----------------------------------------------------------------------
/// In addition to constructing the object using the microwindow
/// ranges, adds a bad sample mask argument.
///
/// The bad sample mask is sized num_bands x num_samples. A value of
/// true in the array means the sample is marked as bad.
//-----------------------------------------------------------------------

template<class T>
SpectralWindowRange::SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges,
        const Array<T, 2>& Bad_sample_mask)
    : range_(Microwindow_ranges)
{
  if(range_.value.depth() != 2) {
    Exception e;
    e << "Microwindow_ranges needs to have a depth of 2\n"
      << "(for lower and upper bound). Found depth of " << range_.value.depth();
    throw e;
  }
  
  // Initialize boolean object from any type that can be evaluated for truth
  bad_sample_mask_.resize(Bad_sample_mask.shape());
  
  for(int i = 0; i < bad_sample_mask_.rows(); ++i) {
    bad_sample_mask_(i, Range::all()) = where(Bad_sample_mask(i, Range::all()), true, false);
  }
}

// See base class for a description.
SpectralBound SpectralWindowRange::spectral_bound() const
{
  Range ra = Range::all();
  std::vector<DoubleWithUnit> lbound, ubound;
  for(int i = 0; i < number_spectrometer(); ++i) {
    DoubleWithUnit lv(min(range_.value(i, ra, ra)), range_.units);
    DoubleWithUnit uv(max(range_.value(i, ra, ra)), range_.units);
    if(range_.units.is_commensurate(units::sample_index) && disp_.size() > 0) {
      SpectralDomain pgrid = disp_[i]->pixel_grid();
      // Special handling for SampleIndex, so it isn't out of range
      if(lv.units.is_commensurate(units::sample_index)) {
	if(lv.value < 1)
	  lv.value = 1;
	if(lv.value > pgrid.data().rows())
	  lv.value = pgrid.data().rows();
	if(uv.value < 1)
	  uv.value = 1;
	if(uv.value > pgrid.data().rows())
	  uv.value = pgrid.data().rows();
      }
      lv = lv.convert_wave(pgrid.units(), pgrid);
      uv = uv.convert_wave(pgrid.units(), pgrid);
      // Might switch the direction of upper and lower
      if(lv.value > uv.value) 
          std::swap(lv.value, uv.value);
    }
    lbound.push_back(lv);
    ubound.push_back(uv);
  }
  return SpectralBound(lbound, ubound);
}

// See base class for a description.

std::vector<int>
SpectralWindowRange::grid_indexes(const SpectralDomain& Grid, int Spec_index) const
{
  range_check(Spec_index, 0, number_spectrometer());
  std::vector<int> res;
  Array<double, 1> g;

  if(range_.units.is_commensurate(units::sample_index)) {
    g.resize(Grid.sample_index().shape());
    g = blitz::cast<double>(Grid.sample_index());
  } else if(range_.units.is_commensurate(units::inv_cm))
    g.reference(Grid.wavenumber(range_.units));
  else
    g.reference(Grid.wavelength(range_.units));
  
  for(int i = 0; i < g.rows(); ++i) {
    bool ok = false;
    
    // Check if sample is bad first, if it is not then check if
    // it falls within the spectral ranges
    if(bad_sample_mask_.rows() == 0 || !bad_sample_mask_(Spec_index, i)) {
      for(int j = 0; j < range_.value.cols(); ++j) {
        if(g(i) >= range_.value(Spec_index, j, 0) &&
           g(i) < range_.value(Spec_index, j, 1)) {
          ok = true;
        }
      }
    }
    
    if(ok) {
      res.push_back(i);
    }
  }

  return res;
}

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void SpectralWindowRange::print(std::ostream& Os) const
{
  Os << "SpectralWindowRange:\n"
     << "   Units: " << range_.units;
  
  for(int i = 0; i < range_.value.rows(); ++i) {
    Os << "  Spec[" << i << "]:\n";
    
    for(int j = 0; j < range_.value.cols(); ++j)
      Os << "    Microwindow[" << j << "]: ("
	 << range_.value(i, j, 0) << ", " << range_.value(i, j, 1) << ")\n";
  }
}
