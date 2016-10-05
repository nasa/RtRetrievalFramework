#ifndef SPECTRAL_DOMAIN_H
#define SPECTRAL_DOMAIN_H
#include "printable.h"
#include "unit.h"
#include "array_ad.h"
#include "array_with_unit.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  For different instruments, it is more natural to either work with
  wavenumbers (e.g., GOSAT) or wavelength (e.g., OCO). Most of our
  code doesn't care if we are using wavenumber or wavelength, so we
  have this one class that can be either. For code that needs one or
  the other, we supply conversion routines to present the data as
  either wavenumber or wavelength.

  As far as I can determine, there isn't any commonly used name that
  means "either wavelength or wavenumber". We've named this
  "SpectralDomain", where "Domain" is used like "Domain and Range" of
  a function, i.e., this is the X-axis of a spectral plot. Perhaps a
  better name will arise and we can rename this class.

  This class is essentially just a blitz::array with units, and the
  additional functionality to convert to wavenumber or wavelength.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/
class SpectralDomain: public Printable<SpectralDomain> {
public:
  enum TypePreference {PREFER_WAVENUMBER, PREFER_WAVELENGTH};

  SpectralDomain(const ArrayAd<double, 1>& Data, 
		 const Unit& U);
  SpectralDomain(const blitz::Array<double, 1>& Data,
		 const Unit& Units = units::inv_cm);
  SpectralDomain(const ArrayWithUnit<double, 1>& Data);
  SpectralDomain(const ArrayAd<double, 1>& Data, 
		 const blitz::Array<int, 1>& Sindex,
		 const Unit& U);
  SpectralDomain(const blitz::Array<double, 1>& Data,
		 const blitz::Array<int, 1>& Sindex,
		 const Unit& Units = units::inv_cm);
  SpectralDomain(const ArrayWithUnit<double, 1>& Data,
		 const blitz::Array<int, 1>& Sindex);

//-----------------------------------------------------------------------
/// Return data. This is either wavenumber or wavelength. This member
/// function is intended for classes that don't care which one we are
/// using. If you do care, then you should call either wavenumber or
/// wavelength. 
///
/// Note that this is a reference to the actual data, so if you intend
/// on modifying this you should make a deep copy.
//-----------------------------------------------------------------------

  const blitz::Array<double, 1>& data() const { return data_.value();}

//-----------------------------------------------------------------------
/// Underlying data, possibly with a Jacobian. The jacobian may have
/// size 0. 
//-----------------------------------------------------------------------

  const ArrayAd<double, 1>& data_ad() const {return data_; }

//-----------------------------------------------------------------------
/// Return sample index. This may be empty if we aren't dealing with
/// the lower resolution grid that maps to a sample index. If present,
/// this gives the sample index for each of the data() points.
///
/// Note that by convention the sample index is 1 based, so the first 
/// sample index is 1 rather than zero.
//-----------------------------------------------------------------------

  const blitz::Array<int, 1>& sample_index() const { return sindex_; }

//-----------------------------------------------------------------------
/// Units that go with data()
//-----------------------------------------------------------------------

  const Unit units() const {return units_;}

//-----------------------------------------------------------------------
/// Clones object into a new copy
//-----------------------------------------------------------------------

  const SpectralDomain clone() const { return SpectralDomain(data_.copy(), units_); }

//-----------------------------------------------------------------------
/// Indicate if this class prefers wavelength or wavenumber. This is
/// what data() is.
//-----------------------------------------------------------------------

  TypePreference type_preference() const
  {
    return (units_.is_commensurate(units::inv_cm) ? 
	    PREFER_WAVENUMBER : PREFER_WAVELENGTH); 
  }

  blitz::Array<double, 1> convert_wave(const Unit& Units) const;
  blitz::Array<double, 1> wavenumber(const Unit& Units = units::inv_cm) const;
  blitz::Array<double, 1> wavelength(const Unit& Units = units::micron) const;
  ArrayWithUnit<double, 1> photon_to_radiance_factor() const;
  void print(std::ostream& Os) const { Os << "SpectralDomain";}

  /// Default constructor needed for SWIG
  SpectralDomain() {}

private:
  ArrayAd<double, 1> data_;
  blitz::Array<int, 1> sindex_;
  Unit units_;
};
}
#endif
