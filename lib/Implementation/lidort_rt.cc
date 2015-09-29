#include "lidort_rt.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(LidortRt, RadiativeTransfer)
.def(luabind::constructor<const boost::shared_ptr<RtAtmosphere>&,
                          const boost::shared_ptr<StokesCoefficient>&,
                          const blitz::Array<double, 1>&,
                          const blitz::Array<double, 1>&, 
                          const blitz::Array<double, 1>&, 
                          bool, int, int, bool>())
.enum_("constants")
[
 luabind::value("maxmoments_input", Lidort_Pars::instance().maxmoments_input)
]
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Atm The atmpsphere to use
/// \param Stokes_coef The Stokes coefficients.
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
/// \param Pure_nadir Flag for controlling azimuth dependence in the
///      output LIDORT will complain if user zenith is 0 and this is
///      not set, when not using ss correction mode
/// \param Number_streams Number of streams to use
/// \param Number_moments Number of moments to use in phase function (-1 to use
///     full resolution
/// \param Do_multi_scatt_only Whether to do multiple scattering calculations
///     only
//-----------------------------------------------------------------------

LidortRt::LidortRt(const boost::shared_ptr<RtAtmosphere>& Atm,
                   const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
                   const blitz::Array<double, 1>& Sza,
                   const blitz::Array<double, 1>& Zen,
                   const blitz::Array<double, 1>& Azm,
                   bool Pure_nadir,
                   int Number_streams,
                   int Number_moments,
                   bool Do_multi_scatt_only)
: SpurrRt(Atm, Stokes_coef, Sza, Zen, Azm)
{   
  rt_driver_.reset(new LidortRtDriver(Number_streams, Number_moments, Do_multi_scatt_only, surface_type(), Zen, Pure_nadir));
}

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void LidortRt::print(std::ostream& Os, bool Short_form) const
{
  Os << "LidortRt\n";
  OstreamPad opad1(Os, "  ");
  SpurrRt::print(opad1, Short_form);
  opad1.strict_sync();
  OstreamPad opad(Os, "    ");
  Os << "  Number stream: " << number_stream() << "\n"
     << "  Number moments: " << number_moment() << "\n"
     << "  Use first order scattering calc: " 
     << (rt_driver()->do_multi_scatt_only() ? "True\n" : "False\n")
     << (rt_driver()->pure_nadir() ? "True\n" : "False\n");

  opad << "\n";
  opad.strict_sync();
}

