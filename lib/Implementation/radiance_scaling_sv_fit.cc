#include "radiance_scaling_sv_fit.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(RadianceScalingSvFit, InstrumentCorrection)
.def(luabind::constructor<const blitz::Array<double, 1>&, 
                          const blitz::Array<bool, 1>&,
                          const DoubleWithUnit&,
                          const std::string&>())
REGISTER_LUA_END()
#endif

void RadianceScalingSvFit::apply_correction
(const SpectralDomain& Pixel_grid,
 const std::vector<int>& Pixel_list,
 SpectralRange& Radiance) const
{
  ArrayAd<double, 1> grid_ad(Pixel_list.size(), Pixel_grid.data_ad().number_variable());
  for (int i = 0; i < (int) Pixel_list.size(); i++) {
    grid_ad(i) = Pixel_grid.data_ad()(Pixel_list[i]);
  }
  grid_ad.resize_number_variable(coeff.number_variable());
  SpectralDomain grid_sd(grid_ad, Pixel_grid.units());
  apply_scaling(grid_sd, Radiance);
}

boost::shared_ptr<InstrumentCorrection> RadianceScalingSvFit::clone() const
{
  return boost::shared_ptr<InstrumentCorrection>
    (new RadianceScalingSvFit(coeff.value().copy(), used_flag, band_ref, band_name));
}

std::string RadianceScalingSvFit::state_vector_name_i(int i) const
{ 
  std::string res ="Radiance Scaling " + band_name;
  if(i == 0)
    res += " Scale";
  else if(i == 1)
    res += " Slope, Ref: " + boost::lexical_cast<std::string>(band_ref.value) + " " + band_ref.units.name();
  else
    res += " Coefficient " + boost::lexical_cast<std::string>(i + 1);
  return res;  
}
 
void RadianceScalingSvFit::print(std::ostream& Os) const
{
  Os << "RadianceScalingSvFit:" << std::endl;
  OstreamPad opad(Os, "    ");
  RadianceScaling::print(opad);
  opad.strict_sync();
}
