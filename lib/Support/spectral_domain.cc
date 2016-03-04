#include "spectral_domain.h"
#include "fp_exception.h"
#include "old_constant.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

std::string spec_domain_type_preference(const SpectralDomain& dom)
{
  if (dom.type_preference() == SpectralDomain::PREFER_WAVENUMBER)
    return "wavenumber";
  else if (dom.type_preference() == SpectralDomain::PREFER_WAVELENGTH)
    return "wavelength";
  else
    return "UNKNOWN";
}

REGISTER_LUA_CLASS(SpectralDomain)
.def("data", &SpectralDomain::data)
.def("units", &SpectralDomain::units)
.def("type_preference", &spec_domain_type_preference)
.def(luabind::constructor<const ArrayAd<double, 1>&, const Unit&>())
.def(luabind::constructor<const Array<double, 1>&>())
.def(luabind::constructor<const Array<double, 1>&, const Unit&>())
.def(luabind::constructor<const ArrayWithUnit<double, 1>&>())
.def(luabind::constructor<const ArrayAd<double, 1>&, const blitz::Array<int, 1>&,const Unit&>())
.def(luabind::constructor<const Array<double, 1>&, const blitz::Array<int, 1>&>())
.def(luabind::constructor<const Array<double, 1>&, const blitz::Array<int, 1>&, const Unit&>())
.def(luabind::constructor<const ArrayWithUnit<double, 1>&, const blitz::Array<int, 1>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
SpectralDomain::SpectralDomain(const ArrayAd<double, 1>& Data,
			       const Unit& Units)
: data_(Data),
  units_(Units)
{
  if(!Units.is_commensurate(units::micron) &&
     !Units.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << Units;
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
SpectralDomain::SpectralDomain(const ArrayAd<double, 1>& Data,
			       const blitz::Array<int, 1>& Sindex,
			       const Unit& Units)
: data_(Data),
  sindex_(Sindex),
  units_(Units)
{
  if(!Units.is_commensurate(units::micron) &&
     !Units.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << Units;
    throw e;
  }
  if(data_.rows() != sindex_.rows())
    throw Exception("Data and Sindex must be the same size");
}

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
 SpectralDomain::SpectralDomain(const blitz::Array<double, 1>& Data,
			       const Unit& Units)
: data_(Data),
  units_(Units)
{
  if(!Units.is_commensurate(units::micron) &&
     !Units.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << Units;
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
 SpectralDomain::SpectralDomain(const blitz::Array<double, 1>& Data,
				const blitz::Array<int, 1>& Sindex,
				const Unit& Units)
: data_(Data),
  sindex_(Sindex),
  units_(Units)
{
  if(!Units.is_commensurate(units::micron) &&
     !Units.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << Units;
    throw e;
  }
  if(data_.rows() != sindex_.rows())
    throw Exception("Data and Sindex must be the same size");
}

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
SpectralDomain::SpectralDomain(const ArrayWithUnit<double, 1>& Data)
  : data_(Data.value),
    units_(Data.units)
{
  if(!units_.is_commensurate(units::micron) &&
     !units_.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << units_;
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Constructor. The units indicate if the data passed in is
/// wavenumber (1/Length) or wavelength (Length)
//-----------------------------------------------------------------------
SpectralDomain::SpectralDomain(const ArrayWithUnit<double, 1>& Data,
			       const blitz::Array<int, 1>& Sindex)
  : data_(Data.value),
    sindex_(Sindex),
    units_(Data.units)
{
  if(!units_.is_commensurate(units::micron) &&
     !units_.is_commensurate(units::inv_cm)) {
    Exception e;
    e << "The units passed to SpectralDomain must be Length or Length^-1.\n"
      << "Units passed in: \n" << units_;
    throw e;
  }
  if(data_.rows() != sindex_.rows())
    throw Exception("Data and Sindex must be the same size");
}

//-----------------------------------------------------------------------
/// Return data as the supplied the units.
//-----------------------------------------------------------------------

Array<double, 1> SpectralDomain::convert_wave(const Unit& Units) const
{
  if(units_.is_commensurate(Units))
    return Array<double, 1>(FullPhysics::conversion(units_, Units) * data_.value());
  else
    return Array<double, 1>(FullPhysics::conversion(1 / units_, Units) / data_.value());
}

//-----------------------------------------------------------------------
/// Return data as wavenumbers. You can optionally supply the units to
/// use.
/// Throws an error if the the optionally supplied units are not
/// commensurate with cm^-1
//-----------------------------------------------------------------------

Array<double, 1> SpectralDomain::wavenumber(const Unit& Units) const
{
  if (Units.is_commensurate(units::inv_cm))
    return convert_wave(Units);
  else {
    stringstream err_msg;
    err_msg << "Supplied units: " 
	    << Units.name() 
	    << " are not commensurate with target units: " 
	    << units::inv_cm.name();
    throw Exception(err_msg.str());
  }
}

//-----------------------------------------------------------------------
/// Return data as wavelengths You can optionally supply the units to
/// use. 
/// Throws an error if the the optionally supplied units are not
/// commensurate with microns
//-----------------------------------------------------------------------

Array<double, 1> SpectralDomain::wavelength(const Unit& Units) const
{
  if (Units.is_commensurate(units::micron))
    return convert_wave(Units);
  else {
    stringstream err_msg;
    err_msg << "Supplied units: " 
	    << Units.name() 
	    << " are not commensurate with target units: " 
	    << units::micron.name();
    throw Exception(err_msg.str());
  } 
}

//-----------------------------------------------------------------------
/// We may want to convert from photon number per second to radiance
/// units. This gives the factor to use in converting.
//-----------------------------------------------------------------------

ArrayWithUnit<double, 1> SpectralDomain::photon_to_radiance_factor() const
{
  using namespace units;
  // Desired output units. Note this is just a convenience, since we
  // typically give wavelength in microns and wavenumber in cm^-1. The
  // conversion could also be handled outside of this class.
  const Unit output_units = (W / inv_cm) / (ph / s / micron);
  const Unit wavenumber_unit = inv_cm;
  const DoubleWithUnit alpha = OldConstant::speed_of_light * OldConstant::planck;
  ArrayWithUnit<double, 1> res;
  res.value.reference
    (Array<double, 1>(alpha.value / wavenumber(wavenumber_unit)));
  res.units = alpha.units / wavenumber_unit / ph;
  return res.convert(output_units);
}
