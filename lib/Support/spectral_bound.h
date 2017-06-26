#ifndef SPECTRAL_BOUND_H
#define SPECTRAL_BOUND_H
#include "printable.h"
#include <vector>
#include "fp_exception.h"
#include "double_with_unit.h"
#include "array_with_unit.h"

namespace FullPhysics {
/****************************************************************//**
  This gives the upper and lower bounds of the SpectralWindow.

  Lower_bound and upper_bound are just meant to give a rough idea of
  the band covered by a spectral window, there is no guarantee that
  lower_bound() + delta or upper_band() - delta will pass the
  grid_indexes test of the SpectralWindow. But the reverse is
  guaranteed, any value < lower_bound or >= upper_bound will certainly
  not pass grid_indexes test.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/

class SpectralBound : public Printable<SpectralBound> {
public:
//-----------------------------------------------------------------------
/// Default constructor.
//-----------------------------------------------------------------------

  SpectralBound() { }

  SpectralBound(const std::vector<DoubleWithUnit>& Lower_bound,
		const std::vector<DoubleWithUnit>& Upper_bound);
  SpectralBound(const ArrayWithUnit<double, 2>& Bound);
  virtual ~SpectralBound() {}

//-----------------------------------------------------------------------
/// Number of spectrometers.
//-----------------------------------------------------------------------

  int number_spectrometer() const { return (int) lower_b.size(); }

//-----------------------------------------------------------------------
/// Center between lower_bound and upper_bound. Turns out we need this
/// often enough to be worth wrapping in a function.
//-----------------------------------------------------------------------

  DoubleWithUnit center(int Spec_index, const Unit& U) const
  {
    return ((lower_bound(Spec_index) + upper_bound(Spec_index)) / 2.0).
      convert_wave(U);
  }

//-----------------------------------------------------------------------
/// Lower bound of window slot.
//-----------------------------------------------------------------------

  DoubleWithUnit lower_bound(int Spec_index) const 
  {
    range_check(Spec_index, 0, number_spectrometer());
    return lower_b[Spec_index];
  }

//-----------------------------------------------------------------------
/// Lower bound but with a unit conversion first in case the conversion
/// reverses ordering.
//-----------------------------------------------------------------------

  DoubleWithUnit lower_bound(int Spec_index, const Unit& U) const
  {
    range_check(Spec_index, 0, number_spectrometer());
    return DoubleWithUnit(std::min(upper_b[Spec_index].convert_wave(U).value,
				   lower_b[Spec_index].convert_wave(U).value),
			  U);
  }

//-----------------------------------------------------------------------
/// Upper bound of window slot.
//-----------------------------------------------------------------------

  DoubleWithUnit upper_bound(int Spec_index) const 
  {
    range_check(Spec_index, 0, number_spectrometer());
    return upper_b[Spec_index];
  }

//-----------------------------------------------------------------------
/// Upper bound but with a unit conversion first in case the conversion
/// reverses ordering.
//-----------------------------------------------------------------------

  DoubleWithUnit upper_bound(int Spec_index, const Unit& U) const
  {
    range_check(Spec_index, 0, number_spectrometer());
    return DoubleWithUnit(std::max(upper_b[Spec_index].convert_wave(U).value,
				   lower_b[Spec_index].convert_wave(U).value),
			  U);
  }

  int spectral_index(const DoubleWithUnit& W) const;
  virtual void print(std::ostream& Os) const;
private:
  std::vector<DoubleWithUnit> lower_b, upper_b;
};
}
#endif
