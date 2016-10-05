#include "spectral_bound.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
typedef DoubleWithUnit (SpectralBound::*bfunc)(int) const;
REGISTER_LUA_CLASS(SpectralBound)
.def("number_spectrometer", &SpectralBound::number_spectrometer)
.def("lower_bound", ((bfunc) &SpectralBound::lower_bound))
.def("upper_bound", ((bfunc) &SpectralBound::upper_bound))
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

SpectralBound::SpectralBound
(const std::vector<DoubleWithUnit>& Lower_bound,
 const std::vector<DoubleWithUnit>& Upper_bound
 )
:lower_b(Lower_bound), upper_b(Upper_bound)
{
  if(lower_b.size() != upper_b.size())
    throw Exception("Lower bound and upper bound need to be the same size");
  for(std::vector<DoubleWithUnit>::size_type i = 0;
      i < lower_b.size(); ++i)
    if(lower_b[i].value > upper_b[i].convert_wave(lower_b[i].units).value)
      throw Exception("Lower bound needs to be less than upper bound");
}

//-----------------------------------------------------------------------
/// Variation of constructor that takes data as an ArrayWithUnit.
//-----------------------------------------------------------------------

SpectralBound::SpectralBound(const ArrayWithUnit<double, 2>& Bound)
{
  if(Bound.cols() != 2) 
    throw Exception("Bound must have cols of 2");
  for(int i = 0; i < Bound.rows(); ++i) {
    if(Bound(i, 0).value > Bound(i, 1).value)
      throw Exception("Lower bound needs to be less than upper bound");
    lower_b.push_back(Bound(i, 0));
    upper_b.push_back(Bound(i, 1));
  }
}

//-----------------------------------------------------------------------
/// Determine spectral index for given wavenumber/wavelength, or
/// return -1 if it doesn't fit in any spectral index.
//-----------------------------------------------------------------------

int SpectralBound::spectral_index(const DoubleWithUnit& W) const
{
  for(int i = 0; i < number_spectrometer(); ++i)
    if(W.value >= lower_bound(i, W.units).value &&
       W.value < upper_bound(i, W.units).value)
      return i;
  return -1;
}

void SpectralBound::print(std::ostream& Os) const
{
  Os << "SpectralBound:\n"
     << "  Number spectrometer: " << number_spectrometer() << "\n";
  for(int i = 0; i < number_spectrometer(); ++i)
    Os << "    [" << i << "]: (" << lower_b[i] << ", "
       << upper_b[i] << ")\n";
}
