#ifndef SPECTRAL_RANGE_H
#define SPECTRAL_RANGE_H
#include "printable.h"
#include "unit.h"
#include "array_with_unit.h"
#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
  We have a number of different spectrums that appear in different
  parts of the code. The spectrum may represent radiances, solar
  spectrum, solar absorption spectrum, etc. In addition to different
  units, the value may have a Jacobian associated with it (e.g., the
  results for RadiativeTransfer), or an uncertainty (e.g., the Level
  1b data). For many purposes, it is convenient to treat these as
  essentially the same thing.

  This class captures this behavior. The data can be accessed just as
  data, as data possibly with a Jacobian, and possibly with an
  associated uncertainty. The data will always be present, but
  depending on the type of Spectrum the uncertainty or jacobian may be
  zero size arrays, indicating they aren't present.

  Similar to SpectralDomain, there doesn't seem to be a commonly used
  name for "stuff that is the Y-axis of a spectral plot". We use the
  name "SpectralRange" where "Range" is used like "Domain and Range"
  of a function. Perhaps a better name will arise and we can rename
  this class.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/
class SpectralRange: public Printable<SpectralRange> {
public:
  SpectralRange(const ArrayWithUnit<double, 1>& Data)
    : data_(Data.value), units_(Data.units) {}
  SpectralRange(const ArrayAd<double, 1>& Data, 
		const Unit& U)
    : data_(Data), units_(U) {}
  SpectralRange(const blitz::Array<double, 1>& Data, 
		const Unit& U)
    : data_(Data), units_(U) {}
  SpectralRange(const blitz::Array<double, 1>& Data, 
		const Unit& U,
		const blitz::Array<double, 1>& Uncertainty)
    : data_(Data), units_(U), 
      uncertainty_(Uncertainty) {}
  SpectralRange(const ArrayAd<double, 1>& Data, 
		const Unit& U,
		const blitz::Array<double, 1>& Uncertainty)
    : data_(Data), units_(U), 
      uncertainty_(Uncertainty) {}

//-----------------------------------------------------------------------
/// Assignment operator so internals are correctly set 
//-----------------------------------------------------------------------
  
  SpectralRange& operator=(const SpectralRange& V)
  { 
    data_.reference(V.data_ad().copy());
    units_ = V.units();
    uncertainty_.reference(V.uncertainty().copy());
    return *this;
  }

//-----------------------------------------------------------------------
/// Underlying data.
//-----------------------------------------------------------------------

  const blitz::Array<double, 1>& data() const {return data_.value(); }

//-----------------------------------------------------------------------
/// Underlying data.
//-----------------------------------------------------------------------

  blitz::Array<double, 1>& data() {return data_.value(); }

//-----------------------------------------------------------------------
/// Underlying data, possibly with a Jacobian. The jacobian may have
/// size 0. 
//-----------------------------------------------------------------------

  const ArrayAd<double, 1>& data_ad() const {return data_; }
  

//-----------------------------------------------------------------------
/// Underlying data, possibly with a Jacobian. The jacobian may have
/// size 0. 
//-----------------------------------------------------------------------

  ArrayAd<double, 1>& data_ad() {return data_; }
  
//-----------------------------------------------------------------------
/// Uncertainty. May be size 0 if we don't have an associated
/// uncertainty. 
//-----------------------------------------------------------------------

  const blitz::Array<double, 1>& uncertainty() const {return uncertainty_;}

//-----------------------------------------------------------------------
/// Clones object into a new copy
//-----------------------------------------------------------------------

  const SpectralRange clone() const { return SpectralRange(data_.copy(), units_, uncertainty_.copy()); }

//-----------------------------------------------------------------------
/// Units of data.
//-----------------------------------------------------------------------
  const Unit& units() const {return units_;}
  SpectralRange convert(const Unit& R) const;
  void print(std::ostream& Os) const { Os << "SpectralRange";}

  /// Default constructor needed for SWIG
  SpectralRange() {}

private:
  ArrayAd<double, 1> data_;
  Unit units_;
  blitz::Array<double, 1> uncertainty_;
};
}
#endif
