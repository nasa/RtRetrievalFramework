#ifndef SPECTRUM_SAMPLING_H
#define SPECTRUM_SAMPLING_H
#include "printable.h"
#include "spectral_domain.h"
#include "spectrum.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This determines the sampling of the spectrum that should be used for
  each of the spectrum indexes.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/
class SpectrumSampling : public Printable<SpectrumSampling> {
public:
  virtual ~SpectrumSampling() {}

//-----------------------------------------------------------------------
/// Number of spectrometers we have.
//-----------------------------------------------------------------------

  int number_spectrometer() const { return nspectrometer; }

//-----------------------------------------------------------------------
/// Wavenumbers/Wavelengths to use for the given spectrometer. We pass
/// in the low resolution grid that we are going to generate after the
/// ILS convolution, along with the ILS half width so we can generate
/// the high resolution points needed to supply the ILS.
//-----------------------------------------------------------------------

  virtual SpectralDomain spectral_domain(int spec_index, 
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Ils_half_width) const = 0;

//-----------------------------------------------------------------------
/// The interpolated spectral domain. The default is that this is just
/// the same as spectral_domain, but derived classes can supply a
/// different implementation if it is doing nonuniform sampling.
//-----------------------------------------------------------------------

  virtual SpectralDomain spectral_domain_interpolated(int Spec_index, 
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Ils_half_width) const
  { return spectral_domain(Spec_index, Lowres_grid, Ils_half_width); }

//-----------------------------------------------------------------------
/// Indicate if spectral_domain and spectral_domain_interpolated are
/// different at all.
//-----------------------------------------------------------------------

  virtual bool need_interpolation(int Spec_index) const { return false; }

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const {Os << "SpectrumSampling";}

protected:
//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
  SpectrumSampling(int num_spectrometer) : nspectrometer(num_spectrometer) {}

//-----------------------------------------------------------------------
/// Default constructor, derived classes should set nspectrometer.
//-----------------------------------------------------------------------
  SpectrumSampling() {}

  int nspectrometer;		//< Derived classes should set this.
};
}
#endif
