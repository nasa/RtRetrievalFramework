#include "solar_doppler_shift_l1b.h"
#include "old_constant.h"
using namespace FullPhysics;
using namespace blitz;
#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(SolarDopplerShiftL1b, SolarDopplerShift)
.def(luabind::constructor<const DoubleWithUnit&, const DoubleWithUnit&,bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create a SolarDopplerShiftL1b.
///
/// \param Solar_distance The distance between observed surface and sun
/// \param Solar_relative_velocity The velocity of the sun along the
///    observed sounding location/sun vector. This includes the
///    rotation of the earth. Positive indicates the sun is moving
///    away from the point.
/// \param Apply_doppler_shift Indicates if we should apply the
///    doppler shift or not. If this is false, then we just 
///    pass the wave numbers through unchanged in doppler_stretch. 
//-----------------------------------------------------------------------

SolarDopplerShiftL1b::SolarDopplerShiftL1b(
const DoubleWithUnit& Solar_distance, 
const DoubleWithUnit& Solar_relative_velocity,
bool Apply_doppler_shift
)
: solar_distance_(Solar_distance),
  solar_relative_velocity_(Solar_relative_velocity),
  apply_doppler_shift_(Apply_doppler_shift)
{
  doppler_shift_ = (solar_relative_velocity_ / 
                    OldConstant::speed_of_light).convert(units::dimensionless).value;
}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void SolarDopplerShiftL1b::print(std::ostream& Os) const
{
  Os << "SolarDopplerShiftL1b\n"
     << "  Solar distance: " << solar_distance_.convert(OldConstant::AU) << "\n"
     << "  Doppler shift:  " << doppler_shift_ << "\n"
     << "  Apply Doppler:  " << (apply_doppler_shift_ ? "true" : "false");
}

// See base class for description
SpectralDomain SolarDopplerShiftL1b::doppler_stretch(
const SpectralDomain& Spec_domain) const
{
  if(apply_doppler_shift_) {
    // The correction is either a multiplication for wavenumbers,
    // or a division for wavelength.
    if(Spec_domain.type_preference() == SpectralDomain::PREFER_WAVENUMBER)
      return SpectralDomain(Array<double, 1>(Spec_domain.data() * 
                                             (1 + doppler_shift_)),
                            Spec_domain.units());
    else
      return SpectralDomain(Array<double, 1>(Spec_domain.data() / 
                                             (1 + doppler_shift_)),
                            Spec_domain.units());
  } else
    return Spec_domain;
}



