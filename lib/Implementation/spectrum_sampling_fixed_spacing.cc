#include "spectrum_sampling_fixed_spacing.h"
#include <fstream>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(SpectrumSamplingFixedSpacing, SpectrumSampling)
.def(luabind::constructor<const ArrayWithUnit<double, 1>&>())
REGISTER_LUA_END()
#endif

// See base class for description.
SpectralDomain SpectrumSamplingFixedSpacing::spectral_domain
(int spec_index,
 const SpectralDomain& Lowres_grid, 
 const DoubleWithUnit& Ils_half_width) const
{
  range_check(spec_index, 0, spec_spacing.rows());

  // Units we input and want to output and the units we
  // use for computing the grid
  Unit u_orig = Lowres_grid.units();
  Unit u_comp = spec_spacing.units;

  // Determines how we compute the gridding, use a shorter name
  DoubleWithUnit sp = spec_spacing(spec_index);

  // Make sure low res grid and ils half width have commensurate units
  if (!u_orig.is_commensurate(Ils_half_width.units)) {
    std::stringstream err_msg;
    err_msg << "Low res grid units:" << std::endl
          << u_orig << std::endl
          << "are not commensurate with ils half width units:" << std::endl
          << Ils_half_width.units;
    throw Exception(err_msg.str());
  }

  // Convert to same units as spacing for the purpose of gridding
  // Make sure data is sorted in ascending order.
  std::vector<DoubleWithUnit> lres_conv;
  BOOST_FOREACH(double v, Lowres_grid.data())
    lres_conv.push_back(DoubleWithUnit(v, u_orig).convert_wave(u_comp));
  std::sort(lres_conv.begin(), lres_conv.end());

  // Fill in points so that we cover +=Ils_half_width for each value listed in
  // low resolution spectral domain. By convention, we have the points
  // be an exact multiple of sp. This isn't strictly necessary, but it
  // could be used to reduce interpolation between spectral points
  // in forward model calculations.
  //
  // Because ILS half width and spectral spacing could be specified
  // using different units, we need to do this conversion back and
  // forth to make sure things are applied appropriately

  std::vector<double> hres;
  if(lres_conv.size() > 0) {
    // Pick the smaller point of +/- since due to conversion
    // to the original units and units not being commensurate
    // one or the other could be the larger value
    DoubleWithUnit begval = 
      min( (lres_conv[0].convert_wave(u_orig) - Ils_half_width).
         convert_wave(u_comp), 
         (lres_conv[0].convert_wave(u_orig) + Ils_half_width).
         convert_wave(u_comp) );

    DoubleWithUnit fpoint = floor(begval / sp) * sp;   
    hres.push_back(fpoint.convert_wave(u_orig).value);

    // Value for for first iteration of loop
    int smax = (int) round(fpoint.value / sp).value;

    BOOST_FOREACH(DoubleWithUnit v, lres_conv) {
      DoubleWithUnit minusval = (v.convert_wave(u_orig) - Ils_half_width).
      convert_wave(u_comp);
      DoubleWithUnit plusval  = (v.convert_wave(u_orig) + Ils_half_width).
      convert_wave(u_comp);

      // Pick values appropriately accounting for ordering change
      // due to conversion of units back and forth
      DoubleWithUnit preval = min(minusval, plusval);
      DoubleWithUnit postval = max(minusval, plusval);

      int fmin = (int) floor(preval / sp).value;
      int fmax = (int) ceil(postval / sp).value;
      for(int f = std::max(fmin, smax + 1); f < fmax; ++f) {
        // Convert back to original units
        DoubleWithUnit fval(f * sp.value, u_comp);
        hres.push_back(fval.convert_wave(u_orig).value);

        // For next iteration
        smax = (int) round(fval / sp).value;
      }
    }
  }
  // Make sure the computed grid points are
  // in sorted order
  std::sort(hres.begin(), hres.end());

  // Create result object
  Array<double, 1> dv((int) hres.size());
  std::copy(hres.begin(), hres.end(), dv.begin());
  SpectralDomain fixed_grid = SpectralDomain(dv, u_orig);
 
  return fixed_grid;
}
