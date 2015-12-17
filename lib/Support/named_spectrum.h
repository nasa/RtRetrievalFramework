#ifndef NAMED_SPECTRUM_H
#define NAMED_SPECTRUM_H
#include "spectrum.h"

namespace FullPhysics {
/****************************************************************//**
 Adds name and spec index fields to a Spectrum. Useful for sending
 Spectrum files to output files.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/
class NamedSpectrum: public Spectrum {
public:
//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

  NamedSpectrum(const SpectralDomain& Spec_domain, 
                const SpectralRange& Spec_range, const std::string& Name,
                int Index)
    : Spectrum(Spec_domain, Spec_range), name_(Name), index_(Index) {}

  NamedSpectrum(const Spectrum& Spec, const std::string& Name, int Index)
    : Spectrum(Spec.spectral_domain(), Spec.spectral_range()),
      name_(Name), index_(Index) {}

//-----------------------------------------------------------------------
/// Name that makes this a named spectrum
//-----------------------------------------------------------------------

  virtual const std::string& name() const {return name_;}

//-----------------------------------------------------------------------
/// An reference index for the spectrum, ie a spectrometer index 
//-----------------------------------------------------------------------

  virtual int index() const {return index_;}

  void print(std::ostream& Os) 
  { 
    Os << "NamedSpectrum:" << std::endl
       << "  name: " << name_ << std::endl
       << "  index: " << index_;
  }

  /// Default constructor needed for SWIG
  NamedSpectrum() {}

private:
  std::string name_;
  int index_;
};
}
#endif
