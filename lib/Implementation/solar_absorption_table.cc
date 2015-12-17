#include "solar_absorption_table.h"
#include "fp_exception.h"
#include "linear_algebra.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(SolarAbsorptionTable, SolarAbsorptionSpectrum)
.def(luabind::constructor<const HdfFile&,
     const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Read the given file for the solar absorption spectrum.
//-----------------------------------------------------------------------

SolarAbsorptionTable::SolarAbsorptionTable(
const HdfFile& F,
const std::string& Hdf_group)
: hdf_file_name(F.file_name()),
  hdf_group(Hdf_group)
{
  ArrayWithUnit<double, 1> sdom = F.read_field_with_unit<double,1>
    (Hdf_group + "/wavenumber");
  ArrayWithUnit<double, 1> srange = F.read_field_with_unit<double,1>
    (Hdf_group + "/spectrum");
  if(!srange.units.is_commensurate(Unit("dimensionless")))
    throw Exception("Solar absorption spectrum units need to be dimensionless");
  if(sdom.value.rows() != srange.value.rows())
    throw Exception("wavenumber and spectrum need to be the same size");
  domain_unit = sdom.units;
  Array<double, 1> x = to_c_order(sdom.value);
  Array<double, 1> y = to_c_order(srange.value);
  table = LinearInterpolate<double, double>
    (x.dataFirst(), x.dataFirst() + x.rows(), y.dataFirst());
}

// See base class for description.
Spectrum SolarAbsorptionTable::solar_absorption_spectrum(
const SpectralDomain& spec_domain) const
{
  Array<double, 1> wv = spec_domain.convert_wave(domain_unit);

  Array<double, 1> res;
  res.resize(wv.shape());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = table(wv(i));
  return Spectrum(spec_domain, SpectralRange(res, Unit("dimensionless")));
}

void SolarAbsorptionTable::print(std::ostream& Os) const
{ 
  Os << "SolarAbsorptionTable\n"
     << "  Hdf file name: " << hdf_file_name << "\n"
     << "  Hdf group:     " << hdf_group << "\n";
}
