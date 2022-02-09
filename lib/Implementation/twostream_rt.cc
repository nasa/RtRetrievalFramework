#include "twostream_rt.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(TwostreamRt, RadiativeTransfer)
.def(luabind::constructor<const boost::shared_ptr<RtAtmosphere>&,
                          const boost::shared_ptr<StokesCoefficient>&,
                          const blitz::Array<double, 1>&,
                          const blitz::Array<double, 1>&, 
                          const blitz::Array<double, 1>&>())
.def(luabind::constructor<const boost::shared_ptr<RtAtmosphere>&,
                          const boost::shared_ptr<StokesCoefficient>&,
                          const blitz::Array<double, 1>&,
                          const blitz::Array<double, 1>&, 
                          const blitz::Array<double, 1>&,
                          bool>())
REGISTER_LUA_END()
#endif

// For debugging purposes, it can be useful to dump out the input
// used by this class. We'll leave this code in place, in case we
// need it again, but normally this turned off.
const bool dump_data = false;

//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Atm Atmosphere class with optical properties to model
/// \param Stokes_coef Multiplier for each stokes component
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
/// \param do_fullquadrature false only for comparison against LIDORT 
//-----------------------------------------------------------------------

TwostreamRt::TwostreamRt(const boost::shared_ptr<RtAtmosphere>& Atm,
                         const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
                         const blitz::Array<double, 1>& Sza,
                         const blitz::Array<double, 1>& Zen,
                         const blitz::Array<double, 1>& Azm, 
                         bool do_fullquadrature)
: SpurrRt(Atm, Stokes_coef, Sza, Zen, Azm)
{   
  rt_driver_.reset(new TwostreamRtDriver(atm->number_layer(), surface_type(), do_fullquadrature));
  if(dump_data)
    std::cout << "# Nlayer:\n" << atm->number_layer() << "\n"
              << "# Surface type:\n" << surface_type() << "\n"
              << "# do_fullquadrature:\n" << do_fullquadrature << "\n";
}

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void TwostreamRt::print(std::ostream& Os, bool Short_form) const
{
  Os << "TwostreamRt\n";
  OstreamPad opad1(Os, "  ");
  SpurrRt::print(opad1, Short_form);
  opad1.strict_sync();
  Os << "do_full_quadrature = " << (rt_driver()->do_full_quadrature() ? "true" : "false") << "\n";
  opad1.strict_sync();
}
