#include "dispersion_fit.h"
#include "old_constant.h"
#include "polynomial_eval.h"
#include "linear_algebra.h"
#include <algorithm>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(DispersionFit)
.def(luabind::constructor<const boost::shared_ptr<Level1b>&>())
.def("fit", &DispersionFit::fit)
.def("shift", &DispersionFit::shift)
REGISTER_LUA_END()
#endif

// Old hard coded values:
// line position of strong, isolated solar line
//const static double ABAND_SOLAR_LINE = 12985.16325e0;
// this takes into account the offset of the GOSAT band 1 ILS.
//const static double ABAND_ILS_OFFSET = 0.203e0; 
// search range
//const static double ABAND_DISP_SEARCH_RANGE[] = { 12983.0, 12988.0 };
//const static double DISPERSION_OFFSET_SCALING[] = { 1, 0.48, 0.38 };

// For solar doppler correction
const static double CON = 6.73951496e-03; // parameter for calculating geocentric from geodetic latitude
const static double EARTH_ROT_SPD = 2.0*OldConstant::pi/86164.09054e0; // earth angular rotation frequency [1/sec]
const static double RADIUS = 6378.137e0; // Earth equatorial radius, in km

//-----------------------------------------------------------------------
/// Initialize class with data needed to perform fit
//-----------------------------------------------------------------------

DispersionFit::DispersionFit(const boost::shared_ptr<Level1b>& Level1b)
: l1b(Level1b)
{
}

blitz::Array<double, 2> DispersionFit::fit(const blitz::Array<double, 2> disp_initial, 
					   const DoubleWithUnit& aband_solar_line_location,
					   const DoubleWithUnit& aband_solar_line_width,
					   const DoubleWithUnit& aband_search_width,
					   const DoubleWithUnit& aband_ils_offset,
					   const ArrayWithUnit<double, 1>& offset_scaling) const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Range ra = Range::all();

  // Data where solar line is present
  Array<double, 1> aband_data = l1b->radiance(0).data();

  // True location of solar line to fit
  double wn_s0 = aband_solar_line_location.value + aband_ils_offset.value;

  // Calculate the Speed of Earth center towards sun
  double V_cen = 497.2 * sin((l1b->time(0).frac_day_of_year() - 4.1) / 365.25 * 2 * OldConstant::pi); // simple approximate model, good to ~ 10 m/s.

  // Calculate the rotational component of the solar doppler velocity
  double sza_r = l1b->solar_zenith(0).convert(Unit("rad")).value;
  double saz_r = l1b->solar_azimuth(0).convert(Unit("rad")).value;
  double lat_r = l1b->latitude(0).convert(Unit("rad")).value;
  double geo_lat = atan(tan(lat_r) / (1.0 + CON));
  double Rloc = 1000.0 * RADIUS / sqrt(1.0 + CON * std::pow(sin(geo_lat), 2));
  double V_rot = -EARTH_ROT_SPD * Rloc * sin(sza_r) * cos(saz_r - 0.5 * OldConstant::pi) * cos(geo_lat);

  double solar_doppler_velocity = V_cen + V_rot;
  
  // Modify position of strong solar line to include solar doppler shift
  // takes into account the solar doppler shift
  double wn_s = wn_s0 * (1.0e0 - solar_doppler_velocity / OldConstant::speed_of_light.value);

  // Create aband dispersion array, 1 indexed
  // Coeffs are in reverse order from that Poly1d expects
  Array<double, 1> aband_wn(aband_data.rows());
  Poly1d aband_poly(disp_initial(0, ra).reverse(firstDim)); 
  for(int pix = 0; pix < aband_wn.rows(); pix++)
    aband_wn(pix) = aband_poly(pix+1);

  // (2) CONSTRUCT THE X, Y FUNCTION TO FIT
  //
  // Use intersection of the points above and below range threshold
  // Make sure to sort the resulting indexes since set.intersection
  // does not guarantee anything about ordering
  Array<double, 1>::const_iterator srch_beg =
    std::lower_bound(aband_wn.begin(), aband_wn.end(), aband_solar_line_location.value - aband_search_width.value);
  Array<double, 1>::const_iterator srch_end = 
    std::upper_bound(aband_wn.begin(), aband_wn.end(), aband_solar_line_location.value + aband_search_width.value);
  Range srch_range(srch_beg.position()(0), srch_end.position()(0));
  Array<double, 1> srch_wn(aband_wn(srch_range));
  Array<double, 1> srch_data(aband_data(srch_range));

  if (srch_wn.rows() == 0)
    throw Exception("Could not find any points in the aband dispersion which encompass the fitting search window");

  Array<double, 1> x( srch_wn - wn_s );

  double mxmeas = max(srch_data); // maximum value
  double m2 = max(mxmeas - srch_data);
  Array<double, 1> y( (mxmeas - srch_data) / m2 ); // this should look like a gaussian

  double cont_mean = 0;
  int cont_count = 0;
  for(int idx = 0; idx < y.rows(); idx++) {
    if(y(idx) < 0.1) {
      cont_mean = cont_mean + y(idx);
      cont_count++;
    }
  }
  cont_mean = cont_mean / cont_count;

  y = y - cont_mean; // subtract off the "continuum"
  int pos = maxIndex(y)(0); // pos = index of the maximum value of y
  double fg = -x(pos);  // first-guess value of spectral shift

  // (3) PERFORM A SIMPLE FIT ASSUMING A PERFECTLY LINEAR MODEL
  //     MODEL IS A GAUSSIAN WITH SIGMA=0.2 CM^-1
  //     FIT PARAMETERS ARE AMPLITUDE AND CENTER OF THE GAUSSIAN
  Array<double, 2> K(y.rows(), 2);
  double sig2 = aband_solar_line_width.value;
  Array<double, 1> ygauss( exp(-pow2(x + fg) / (2.0 * sig2)) );

  K(ra, 0) = ygauss;
  K(ra, 1) = -ygauss / sig2 * (x + fg);

  // Note that KtK is a 2x2 symmetric matrix and can be inverted analytically if desired.
  // form matrix (Kt K)^{-1} Kt
  Array<double, 2> Kt( K.transpose(secondDim, firstDim) );

  // Kt * K, matrix multiplication using blitz tensor notaion
  Array<double, 2> KtK( sum(Kt(i1, i3) * K(i3, i2), i3) ); 
  Array<double, 2> KtK_m1 = generalized_inverse(KtK);
  Array<double, 1> yminusg(y - ygauss);
  Array<double, 1> Kty(sum(Kt(i1, i2) * yminusg(i2), i2));
  Array<double, 1> delta(sum(KtK_m1(i1, i2) * Kty(i2), i2));

  Array<double, 1> fit(delta.rows());
  fit = 1, fg;
  fit = fit + delta;

  double aband_shift = fit(1);

  if(0) {
    std::cerr << "#sig2 = " << std::endl
	      << sig2 << std::endl
	      << "# x = " << std::endl
	      << x << std::endl
	      << "# y = " << std::endl
	      << y << std::endl
	      << "# fit = " <<std::endl
	      << fit << std::endl
	      << "# offset_scaling = "<<std::endl
	      << offset_scaling.value << std::endl
	      << "# aband_shift = " << std::endl
	      << aband_shift << std::endl;
    // Output to plot stuff like so:
    // pyplot.plot(x, y) # measured data
    // yfit = fit[0]* numpy.exp(-(x+fit[1])**2/ (2.0*sig2))
    // pyplot.plot(x,yfit) # fit modeled data
  }

  // Create dispersion by scaling non aband values according to this method:
  // Where below wn_spc is the spacing, ie second coefficient
  // wco2 offset : delta_2 = 0.48 * ( delta_1 + wn_spc) - wn_spc)
  // sco2_offset : delta_3 = 0.38 * ( delta_1 + wn_spc) - wn_spc)
  Array<double, 2> disp_adj(disp_initial.copy());

  // Store computed shift for each band
  shift_.resize(disp_adj.rows());

  for(int band_idx = 0; band_idx < disp_adj.rows(); band_idx++) {
    shift_(band_idx) = offset_scaling.value(band_idx) * aband_shift;
    disp_adj(band_idx, 0) = disp_adj(band_idx, 0) + shift_(band_idx);
  }

  return disp_adj;
}

void DispersionFit::print(std::ostream& Os) const
{
  Os << "DispersionFit\n";
}
