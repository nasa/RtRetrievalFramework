#ifndef RADIATIVE_TRANSFER_H
#define RADIATIVE_TRANSFER_H
#include "printable.h"
#include "accumulated_timer.h"
#include "array_ad.h"
#include "spectrum.h"
#include <blitz/array.h>
#include <boost/progress.hpp>
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/****************************************************************//**
  This runs a Radiative Transfer code to determine the reflectance for a
  given set of wavelengths.

  We support both vector and scalar calculations. Because of the large
  size of the arrays returned, we often use only a subset of the
  stokes parameters given by number_stokes(). This can be up to 4, in
  which case we return I, Q, U and V (in that order). For Gosat, we
  commonly return 3 parameters: I, Q and U.

  If the Radiative Transfer code is scalar, then you can either set
  the number_stokes() to 1 and return I, or just set the terms other
  than I to 0.
*******************************************************************/

class RadiativeTransfer : public Printable<RadiativeTransfer> {
public:
  virtual ~RadiativeTransfer() {}

//-----------------------------------------------------------------------
/// Number of stokes parameters we will return in stokes and 
/// stokes_and_jacobian.
//-----------------------------------------------------------------------

  virtual int number_stokes() const = 0;

//-----------------------------------------------------------------------
/// Number of spectrometer we have.
//-----------------------------------------------------------------------

  virtual int number_spectrometer() const = 0;

//-----------------------------------------------------------------------
/// Calculate reflectance for the given set of wavenumbers/wavelengths.
///
/// \param Spec_domain List of wavenumber/wavelength to calculate for.
/// \param Spec_index The Spectral index
/// \param Skip_jacobian If true, don't do the Jacobian
/// calculation. Often this is significantly faster to calculate.
/// \return The set of reflectance values.
//-----------------------------------------------------------------------

  virtual Spectrum reflectance
  (const SpectralDomain& Spec_domain, int Spec_index, 
   bool Skip_jacobian = false) const = 0;

//-----------------------------------------------------------------------
/// Calculate stokes vector for the given set of wavenumbers/wavelengths. 
///
/// \param Spec_domain List of wavenumber/wavelength to calculate for.
/// \param Spec_index The Spectral index
/// \return The set of stokes coefficients. This is
///      Spec_domain.data().rows() x number_stokes() in size.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain,
					 int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Calculate stokes vector for the given set of
/// wavenumbers/wavelengths. This also calculates the Jacobian of the
/// stokes with respect to the state vector elements.
///
/// \param Spec_domain List of wavenumber/wavelength to calculate for.
/// \param Spec_index The Spectral index
/// \return The set of stokes coefficients, along with derivatives with
///   respect to the state vector elements. This is
///   Spec_domain.data().rows() x number_stokes() in size.
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 2> stokes_and_jacobian
  (const SpectralDomain& Spec_domain, int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os, bool Short_form = false) const 
  { Os << "RadiativeTransfer";}
protected:
  static AccumulatedTimer timer;
  boost::shared_ptr<boost::progress_display> progress_display(const 
		   blitz::Array<double, 1>& wn) const;
};

}
#endif
