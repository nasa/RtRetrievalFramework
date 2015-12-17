#ifndef SPECTRUM_H
#define SPECTRUM_H
#include "spectral_domain.h"
#include "spectral_range.h"

namespace FullPhysics {
/****************************************************************//**
  This is a full spectrum, which contains a SpectralRange and
  SpectralDomain. The SpectralRange has a value, possibly with
  associated Jacobian and/or uncertainty, and units. The
  SpectralDomain has the wavenumber/wavelength of the spectrum.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/
class Spectrum: public Printable<Spectrum> {
public:
//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

  Spectrum(const SpectralDomain& Spec_domain, 
	   const SpectralRange& Spec_range)
    : spec_domain_(Spec_domain), spec_range_(Spec_range) {}

//-----------------------------------------------------------------------
/// Spectral domain (i.e., wavenumber or wavelength).
//-----------------------------------------------------------------------

  const SpectralDomain& spectral_domain() const {return spec_domain_;}
  SpectralDomain& spectral_domain() {return spec_domain_;}

//-----------------------------------------------------------------------
/// Spectral range (e.g, radiance values)
//-----------------------------------------------------------------------

  const SpectralRange& spectral_range() const {return spec_range_;}
  SpectralRange& spectral_range() {return spec_range_;}

//-----------------------------------------------------------------------
/// Clones object into a new copy
//-----------------------------------------------------------------------

  const Spectrum clone() const { return Spectrum(spec_domain_.clone(), spec_range_.clone()); }

  void print(std::ostream& Os) { Os << "Spectrum";}

  /// Default constructor needed for SWIG
  Spectrum() {}
  // Put in virtual table, again for use by SWIG
  virtual ~Spectrum() {}
private:
  SpectralDomain spec_domain_;
  SpectralRange spec_range_;
};
}
#endif
