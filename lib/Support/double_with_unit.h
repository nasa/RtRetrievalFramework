#ifndef DOUBLE_WITH_UNIT_H
#define DOUBLE_WITH_UNIT_H
#include "printable.h"
#include "unit.h"
#include <cmath>

namespace FullPhysics {
  class SpectralDomain;
/****************************************************************//**
  We frequently have a double with units associated with it. This is a
  simple structure that just keeps these two things together.
*******************************************************************/
class DoubleWithUnit : public Printable<DoubleWithUnit>,
		       boost::field_operators<DoubleWithUnit> {
public:
  DoubleWithUnit() {}
  DoubleWithUnit(double V, const Unit& U)
    : value(V), units(U) {}
  DoubleWithUnit(double V, const std::string& U)
    : value(V), units(U) {}
  DoubleWithUnit(double V)
    : value(V), units(units::dimensionless) {}
  double value;
  Unit units;

//-----------------------------------------------------------------------
/// Basic math operators for class.
//-----------------------------------------------------------------------
  inline  DoubleWithUnit& operator*=(const DoubleWithUnit& D)
  { value *= D.value; units *= D.units; return *this;}
  inline DoubleWithUnit& operator/=(const DoubleWithUnit& D)
  { value /= D.value; units /= D.units; return *this;}
  inline DoubleWithUnit& operator+=(const DoubleWithUnit& D)
  { value += D.value * FullPhysics::conversion(D.units, units); return *this;}
  inline DoubleWithUnit& operator-=(const DoubleWithUnit& D)
  { value -= D.value * FullPhysics::conversion(D.units, units); return *this;}

//-----------------------------------------------------------------------
/// Convert to the given units.
//-----------------------------------------------------------------------

  inline DoubleWithUnit convert(const Unit& R) const
  { return DoubleWithUnit(value * FullPhysics::conversion(units, R), R); }

//-----------------------------------------------------------------------
/// We often need to handle conversion from wavenumber to/from
/// wavelength. This is either a normal conversion of the units before
/// and after match in the power of length (so cm^-1 to m^-1), or do
/// an inversion. Since we do this often enough, it is worth having a
/// function that handles this logic.
//-----------------------------------------------------------------------
    
  inline DoubleWithUnit convert_wave(const Unit& R) const
  {
    if(units.is_commensurate(R))
      return convert(R);
    else
      return (1.0 / *this).convert(R);
  }

  DoubleWithUnit convert_wave(const Unit& R, 
			      const SpectralDomain& Pixel_grid) const;

  void print(std::ostream& Os) const 
  { Os << value << " " << units.name(); }

};

//-----------------------------------------------------------------------
/// \ingroup Miscellaneous
// Compare DoubleWithUnits
///
/// We define <=, >= and > in terms of this operator.
//-----------------------------------------------------------------------
inline bool operator<(const FullPhysics::DoubleWithUnit& A, const FullPhysics::DoubleWithUnit& B)
{ return A.value < B.value; }

}

//-----------------------------------------------------------------------
/// Math functions.
//-----------------------------------------------------------------------
namespace std {
  // Math functions are in std:: namespace.
  inline FullPhysics::DoubleWithUnit floor(const FullPhysics::DoubleWithUnit& x) 
  {
    return FullPhysics::DoubleWithUnit(::floor(x.value), x.units);
  }

  inline FullPhysics::DoubleWithUnit ceil(const FullPhysics::DoubleWithUnit& x) 
  {
    return FullPhysics::DoubleWithUnit(::ceil(x.value), x.units);
  }

  inline FullPhysics::DoubleWithUnit round(const FullPhysics::DoubleWithUnit& x) 
  {
    return FullPhysics::DoubleWithUnit(::round(x.value), x.units);
  }

  inline FullPhysics::DoubleWithUnit min(const FullPhysics::DoubleWithUnit& x, const FullPhysics::DoubleWithUnit& y) 
  {
    return FullPhysics::DoubleWithUnit(std::min(x.value, y.value), x.units);
  }

  inline FullPhysics::DoubleWithUnit max(const FullPhysics::DoubleWithUnit& x, const FullPhysics::DoubleWithUnit& y) 
  {
    return FullPhysics::DoubleWithUnit(std::max(x.value, y.value), x.units);
  }

}

#endif
