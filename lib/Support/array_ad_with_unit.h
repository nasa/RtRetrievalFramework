#ifndef ARRAY_AD_WITH_UNIT_H
#define ARRAY_AD_WITH_UNIT_H
#include "array_ad.h"
#include "auto_derivative_with_unit.h"
#include "unit.h"

namespace FullPhysics {
/****************************************************************//**
  This is a ArrayAd that also has units associated with it.
  This is a simple structure that just keeps these two things
  together. 
*******************************************************************/

template<class T, int D> class ArrayAdWithUnit
{
public:
  ArrayAdWithUnit() {}
  ArrayAdWithUnit(const ArrayAd<T, D>& V, const Unit& U)
    : value(V), units(U) {}
  ArrayAdWithUnit(const ArrayAd<T, D>& V, const std::string& U)
    : value(V), units(U) {}
  ArrayAdWithUnit(const ArrayAd<T, D>& V)
    : value(V), units(units::dimensionless) {}
  ArrayAd<T, D> value;
  Unit units;

//-----------------------------------------------------------------------
/// Convert to the given units.
//-----------------------------------------------------------------------

  inline ArrayAdWithUnit<T, D> convert(const Unit& R) const
  { ArrayAd<T, D> res(value.copy());
    double c = FullPhysics::conversion(units, R);
    res.value() *= c;
    res.jacobian() *= c;
    return ArrayAdWithUnit<T, D>(res, R);
  }
  AutoDerivativeWithUnit<T> operator()(int i1) const
  { return AutoDerivativeWithUnit<T>(value(i1), units); }
  AutoDerivativeWithUnit<T> operator()(int i1, int i2) const
  { return AutoDerivativeWithUnit<T>(value(i1, i2), units); }
  AutoDerivativeWithUnit<T> operator()(int i1, int i2, int i3) const
  { return AutoDerivativeWithUnit<T>(value(i1, i2, i3), units); }
  AutoDerivativeWithUnit<T> operator()(int i1, int i2, int i3, int i4) const
  { return AutoDerivativeWithUnit<T>(value(i1, i2, i3, i4), units); }
  AutoDerivativeWithUnit<T> operator()(int i1, int i2, int i3, int i4, 
				       int i5) const
  { return AutoDerivativeWithUnit<T>(value(i1, i2, i3, i4, i5), units); }
  int rows() const {return value.rows();}
  int cols() const {return value.cols();}
  int depth() const {return value.depth();}
  bool is_constant() const {return value.is_constant();}
  int number_variable() const {return value.number_variable(); }
  void reference(const ArrayAdWithUnit<T, D>& V)
  {
    value.reference(V.value);
    units = V.units;
  }
};
}
#endif
