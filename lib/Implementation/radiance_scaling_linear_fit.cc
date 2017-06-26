#include "radiance_scaling_linear_fit.h"
#include "ostream_pad.h"
#include "linear_algebra.h"

#include <fstream>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(RadianceScalingLinearFit, InstrumentCorrection)
.def(luabind::constructor<const SpectralRange&, 
                          const DoubleWithUnit&,
                          const std::string&>())
.def(luabind::constructor<const SpectralRange&, 
                          const DoubleWithUnit&,
                          const std::string&,
                          const bool>())
REGISTER_LUA_END()
#endif

void RadianceScalingLinearFit::apply_correction
(const SpectralDomain& Pixel_grid,
 const std::vector<int>& Pixel_list,
 SpectralRange& Radiance) const
{
  ArrayAd<double, 1> grid_ad(Pixel_list.size(), Pixel_grid.data_ad().number_variable());
  for (int i = 0; i < (int) Pixel_list.size(); i++) {
    grid_ad(i) = Pixel_grid.data_ad()(Pixel_list[i]);
  }
  SpectralDomain grid_sd(grid_ad, Pixel_grid.units());

  // Convert measured radiances into same units as radiance to scale
  // and only use those pixels being used for the retrieval
  double conv_factor = FullPhysics::conversion(measured_radiance.units(), Radiance.units());
  Array<double, 1> meas_conv_data(Pixel_list.size());
  for(int i = 0; i < (int) Pixel_list.size(); i++)
    meas_conv_data(i) = measured_radiance.data()(Pixel_list[i]) * conv_factor;
  SpectralRange meas_conv(meas_conv_data, Radiance.units());

  // Formulate radiance scaling as a linear equation and solve
  Array<double, 2> problem_mat;
  if (do_offset)
    problem_mat.resize(Pixel_list.size(), 3);
  else
    problem_mat.resize(Pixel_list.size(), 2);

  double band_ref_conv = band_ref.convert_wave(grid_sd.units()).value;
  Array<double, 1> wrefd( grid_sd.data() - band_ref_conv );

  Range all = Range::all();

  // Create problem matrix
  // r0*x0 + r0*x1*wrefd + x2
  // To solve for scaling coefficients (x0, x1)
  // and possibly an offset x2
  problem_mat(all, 0) = Radiance.data();
  problem_mat(all, 1) = Radiance.data() * wrefd;
  if (do_offset)
    problem_mat(all,2) = 1.0;

  // Use a QR factorization based least square solver since
  // our problem is overdetermined
  Array<double, 1> scale_offset = solve_least_squares_qr(problem_mat, meas_conv.data());

  scaling_coeff(all) = scale_offset(Range(0,1));
  if (do_offset)
    offset = scale_offset(2);
 
  // Debugging output
  if(false) {
    std::ofstream debug_out(("radiance_scaling_linear_fit_debug_" + band_name + ".txt").c_str());
    debug_out << "# scale_offset: " << std::endl << scale_offset << std::endl;
    double conv_factor = FullPhysics::conversion(measured_radiance.units(), Radiance.units());
    debug_out << "# inst_corr conv_factor = " << conv_factor << std::endl
      << "# measured_radiance: " << measured_radiance.units().name() << std::endl
      << measured_radiance.data() << std::endl
      << "# ---" << std::endl
      << "# measured_radiance converted = " << meas_conv.units().name() << std::endl
      << meas_conv.data() << std::endl
      << "# ---" << std::endl
      << "# Radiance: " << Radiance.units().name() << std::endl
      << Radiance.data() << std::endl
      << "# ---" << std::endl
      << "# wrefd: " << std::endl
      << wrefd << std::endl;
    debug_out.close();
  }  

  apply_scaling(grid_sd, Radiance);
}

boost::shared_ptr<InstrumentCorrection> RadianceScalingLinearFit::clone() const
{
  return boost::shared_ptr<InstrumentCorrection>
    (new RadianceScalingLinearFit(measured_radiance, band_ref, band_name, do_offset));
}

void RadianceScalingLinearFit::print(std::ostream& Os) const
{
  Os << "RadianceScalingLinearFit:" << std::endl;
  OstreamPad opad(Os, "    ");
  RadianceScaling::print(opad);
  opad << "Use offset: " << (do_offset ? "True" : "False") << std::endl;
  opad.strict_sync();
}
