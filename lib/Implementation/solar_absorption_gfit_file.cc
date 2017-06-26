#include "solar_absorption_gfit_file.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(SolarAbsorptionGfitFile, SolarAbsorptionSpectrum)
.def(luabind::constructor<const std::string&, double>())
REGISTER_LUA_END()
#endif

extern "C" {
  void solar_pts(const int *lunr, const char* filename, const int* filename_len, const double *fzero, const double *grid, const double *frac, double *spts, const int *ncp);
}

//-----------------------------------------------------------------------
/// Read the given line list file, and use for calculating the solar
/// absorption spectrum.
///
/// \param Line_list_file Line list file
/// \param Fraction_solar_diameter Fraction of Solar diameter viewed.
//-----------------------------------------------------------------------

SolarAbsorptionGfitFile::SolarAbsorptionGfitFile(const std::string& Line_list_file,
                                                 double Fraction_solar_diameter)
: line_list_file_(Line_list_file), fraction_solar_diameter_(Fraction_solar_diameter)
{
}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void SolarAbsorptionGfitFile::print(std::ostream& Os) const
{
  Os << "Solar Absorption GFIT File:\n"
     << "  Fraction solar diameter: " << fraction_solar_diameter_;
}

// See base class for description.
Spectrum SolarAbsorptionGfitFile::solar_absorption_spectrum( const SpectralDomain& spec_domain) const
{
  int lunr = 99;
  // Set up inputs to match our grid
  //  V(i) = fzero + i * grid       i = 1,NCP
  double grid = spec_domain.data()(1) - spec_domain.data()(0);
  double fzero = spec_domain.data()(0) - grid;
  int ncp = spec_domain.data().rows();
  Array<double, 1> spts(ncp);
  int file_len = line_list_file_.length();
  solar_pts(&lunr, line_list_file_.c_str(), &file_len, &fzero, &grid, &fraction_solar_diameter_, spts.dataFirst(), &ncp);

  return Spectrum(spec_domain, SpectralRange(spts, units::dimensionless));
}
