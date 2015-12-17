#include "solar_continuum_table.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include <boost/foreach.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(SolarContinuumTable, SolarContinuumSpectrum)
.def(luabind::constructor<const HdfFile&, const std::string&, bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor. There are two conventions used for the reported solar
/// spectrum. For OCO, we give the results as ph/s/m^2/micron, while
/// for GOSAT we convert using the number of photons at each
/// wavenumber to give W/m^2/cm^-1. You can choose the convention used
/// by specifying if we Convert_from_photon or not.
//-----------------------------------------------------------------------

SolarContinuumTable::SolarContinuumTable
(const HdfFile& F, const std::string& Hdf_group,bool Convert_from_photon)
: convert_from_photon(Convert_from_photon),
  hdf_file_name(F.file_name()),
  hdf_group(Hdf_group)
{
  ArrayWithUnit<double, 1> sdom = F.read_field_with_unit<double,1>
    (Hdf_group + "/wavenumber");
  ArrayWithUnit<double, 1> srange = F.read_field_with_unit<double,1>
    (Hdf_group + "/spectrum");
  if(!srange.units.is_commensurate(Unit("ph / (s * m * m * micron)")))
    throw Exception("Solar continuum spectrum units need to be commensurate with ph / (s * m * m * micron)");
  if(sdom.value.rows() != srange.value.rows())
    throw Exception("wavenumber and spectrum need to be the same size");
  domain_unit = sdom.units;
  range_unit = srange.units;
  Array<double, 1> x = to_c_order(sdom.value);
  Array<double, 1> y = to_c_order(srange.value);
  table = LinearInterpolate<double, double>
    (x.dataFirst(), x.dataFirst() + x.rows(), y.dataFirst());
}

// See base class for description.
Spectrum SolarContinuumTable::solar_continuum_spectrum(
  const SpectralDomain& spec_domain) const
{
  Array<double, 1> wv = spec_domain.convert_wave(domain_unit);

  ArrayWithUnit<double, 1> res;
  res.value.resize(wv.shape());
  res.units = range_unit;
  for(int i = 0; i < res.value.rows(); ++i)
    res.value(i) = table(wv(i));

  // Convert from ph/s/m^2/micron to W/m^2/cm^-1 if requested
  if(convert_from_photon)
    res *= spec_domain.photon_to_radiance_factor();
  return Spectrum(spec_domain, SpectralRange(res.value, res.units));
}

void SolarContinuumTable::print(std::ostream& Os) const
{ 
  Os << "SolarContinuumTable\n"
     << "  Hdf file name:       " << hdf_file_name << "\n"
     << "  Hdf group:           " << hdf_group << "\n"
     << "  Convert from photon: " << (convert_from_photon ? "true" : "false")
     << "\n";
}
