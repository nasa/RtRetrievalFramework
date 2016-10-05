#include "solar_absorption_oco_file.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include <sys/stat.h>

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(SolarAbsorptionOcoFile, SolarAbsorptionSpectrum)
.def(luabind::constructor<const HdfFile&,
     const std::string&>())
REGISTER_LUA_END()
#endif

extern "C" {
  void saof_sunspect(const int* nlines, const double* freq,
		     const double* stren, const double* w_wid,
		     const double* d_wid,
		     const int* ncp, const double* wn_array_c, 
		     const double* frac, const int* units, 
		     double* sot_c);
}

//-----------------------------------------------------------------------
/// Read the given line list file, and use for calculating the solar
/// absorption spectrum.
///
/// \param Hdf_static_input HDF input file
/// \param Hdf_group HDF group to get data from
/// \param Fraction_solar_diameter Fraction of Solar diameter viewed.
///
/// \todo Note that the Fraction_solar_diameter parameter is mostly
/// ignored, the underlying Fortran code ignores this value and sets
/// the variable "sld" to 1.0 regardless of the Fraction_solar_diameter
/// passed in.
//-----------------------------------------------------------------------

SolarAbsorptionOcoFile::SolarAbsorptionOcoFile(
     const HdfFile& Hdf_static_input,
     const std::string& Hdf_group,
     double Fraction_solar_diameter)
: hdf_file_name(Hdf_static_input.file_name()),
  hdf_group(Hdf_group),
  fraction_solar_diameter_(Fraction_solar_diameter)
{
  freq.reference(Hdf_static_input.read_field_with_unit<double, 1>
  		 (Hdf_group + "/Solar_Line_List/frequency").
  		 convert(units::inv_cm).value);
  stren.reference(Hdf_static_input.read_field<double, 1>
  		 (Hdf_group + "/Solar_Line_List/optical_thickness"));
  w_wid.reference(Hdf_static_input.read_field_with_unit<double, 1>
  		 (Hdf_group + "/Solar_Line_List/folding_width").
  		  convert(units::inv_cm).value);
  d_wid.reference(Hdf_static_input.read_field_with_unit<double, 1>
  		 (Hdf_group + "/Solar_Line_List/doppler_width").
  		  convert(units::inv_cm).value);
  // stren should not be less that 0.
  stren = where(stren < 0, 0.0, stren);
  // Make sure frequency is sorted. This is required by the Fortran
  // code.
  for(int i = 1; i < freq.rows(); ++i)
    if(freq(i) < freq(i-1)) {
      throw Exception("The frequency needs to be sorted in the solar line list.");
    }
}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void SolarAbsorptionOcoFile::print(std::ostream& Os) const
{
  Os << "Solar Absorption OCO File:\n"
     << "  Hdf file name:           " << hdf_file_name << "\n"
     << "  Hdf group:               " << hdf_group << "\n"
     << "  Fraction solar diameter: " << fraction_solar_diameter_;
}

// See base class for description.
Spectrum SolarAbsorptionOcoFile::solar_absorption_spectrum(
const SpectralDomain& spec_domain) const
{
  // wave_number might not be contiguous. Since the Fortran code
  // requires that it is, we make a copy if needed.
  blitz::Array<double, 1> wn = 
    to_fortran_const(spec_domain.wavenumber(units::inv_cm));
  int number_wavelength = wn.rows();
  blitz::Array<double,1> res(wn.shape());
  res = 0;
  // "0" here refers to the units of the wavelength. We always use
  // wavenumbers, which is "0"
  int units = 0;
  int nline = number_line();
  if(number_wavelength >= 1)	// Empty data is ok, and is a useful
				// edge case. But Fortran can't handle
				// this, so skip call.
    saof_sunspect(&nline, freq.dataFirst(), stren.dataFirst(),
		  w_wid.dataFirst(), d_wid.dataFirst(),
		  &number_wavelength,
		  wn.dataFirst(),
		  &fraction_solar_diameter_, &units, res.dataFirst());
  blitz::Array<double, 1> res_exp(exp(res));
  return Spectrum(spec_domain, 
		  SpectralRange(res_exp, units::dimensionless));
}
