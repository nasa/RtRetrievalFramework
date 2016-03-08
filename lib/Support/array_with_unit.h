#ifndef ARRAY_WITH_UNIT_H
#define ARRAY_WITH_UNIT_H
#include "printable.h"
#include "unit.h"
#include "double_with_unit.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  We frequently have a array of numbers with units associated with
  them. This is a simple structure that just keeps these two things
  together. 
*******************************************************************/
template<class T, int D> 
class ArrayWithUnit : public Printable<ArrayWithUnit<T, D> >,
                             boost::field_operators<ArrayWithUnit<T, D> > {
public:
  ArrayWithUnit() { }
  ArrayWithUnit(const ArrayWithUnit<T, D>& V)
    : value(V.value.copy()), units(V.units) { }
  ArrayWithUnit(const blitz::Array<T, D>& Value, const Unit& Value_units)
    : value(Value.copy()), units(Value_units) { }
  ArrayWithUnit(const blitz::Array<T, D>& Value, const std::string& Value_units_name)
    : value(Value.copy()), units(Value_units_name) { }

  blitz::Array<T, D> value;
  Unit units;
  
//-----------------------------------------------------------------------
/// Assignment operator so internals are correctly set 
//-----------------------------------------------------------------------
  ArrayWithUnit<T, D>& operator=(const ArrayWithUnit<T, D>& V)
  { value.reference(V.value.copy()); units = V.units; return *this;}
  
//-----------------------------------------------------------------------
/// Basic math operators for class.
//-----------------------------------------------------------------------
  inline ArrayWithUnit<T, D>& operator*=(const ArrayWithUnit<T, D>& V)
  { value *= V.value; units *= V.units; return *this;}
  inline ArrayWithUnit<T, D>& operator/=(const ArrayWithUnit<T, D>& V)
  { value /= V.value; units /= V.units; return *this;}
  inline ArrayWithUnit<T, D>& operator+=(const ArrayWithUnit<T, D> & V)
  { value += V.value * FullPhysics::conversion(V.units, units); return *this;}
  inline ArrayWithUnit<T, D>& operator-=(const ArrayWithUnit<T, D> & V)
  { value -= V.value * FullPhysics::conversion(V.units, units); return *this;}

  DoubleWithUnit operator()(int i1) const 
  {return DoubleWithUnit(value(i1), units);}
  DoubleWithUnit operator()(int i1, int i2) const 
  {return DoubleWithUnit(value(i1, i2), units);}
  DoubleWithUnit operator()(int i1, int i2, int i3) const 
  {return DoubleWithUnit(value(i1, i2, i3), units);}
  DoubleWithUnit operator()(int i1, int i2, int i3, int i4) const 
  {return DoubleWithUnit(value(i1, i2, i3, i4), units);}

  // Convenience wrappers so we don't need to enumerate
  // all possible instances where say an int and Range
  // are interchanged
  template<class I>
  ArrayWithUnit<T, D> operator()(I i1) const 
  {return ArrayWithUnit<T, D>(value(i1), units);}

  template<class I>
  ArrayWithUnit<T, D> operator()(I i1, I i2) const 
  {return ArrayWithUnit<T, D>(value(i1, i2), units);}
  template<class I>
  ArrayWithUnit<T, D> operator()(I i1, I i2, I i3) const 
  {return ArrayWithUnit<T, D>(value(i1, i2, i3), units);}
  template<class I>
  ArrayWithUnit<T, D> operator()(I i1, I i2, I i3, I i4) const 
  {return ArrayWithUnit<T, D>(value(i1, i2, i3, i4), units);}

  template<class I, class J>
  ArrayWithUnit<T, D-1> operator()(I i1, J i2) const 
  {return ArrayWithUnit<T, D-1>(value(i1, i2), units);}

  template<class I, class J>
  ArrayWithUnit<T, D-1> operator()(J i1, I i2, I i3) const 
  {return ArrayWithUnit<T, D-1>(value(i1, i2, i3), units);}
  template<class I, class J>
  ArrayWithUnit<T, D-1> operator()(I i1, J i2, I i3) const 
  {return ArrayWithUnit<T, D-1>(value(i1, i2, i3), units);}
  template<class I, class J>
  ArrayWithUnit<T, D-1> operator()(I i1, I i2, J i3) const 
  {return ArrayWithUnit<T, D-1>(value(i1, i2, i3), units);}

  template<class I, class J>
  ArrayWithUnit<T, D-2> operator()(I i1, J i2, J i3) const 
  {return ArrayWithUnit<T, D-1>(value(i1, i2, i3), units);}
  template<class I, class J>
  ArrayWithUnit<T, D-2> operator()(J i1, I i2, J i3) const 
  {return ArrayWithUnit<T, D-1>(value(i1, i2, i3), units);}
  template<class I, class J>
  ArrayWithUnit<T, D-2> operator()(J i1, J i2, I i3) const 
  {return ArrayWithUnit<T, D-1>(value(i1, i2, i3), units);}

//-----------------------------------------------------------------------
/// Convert to the given units. 
//-----------------------------------------------------------------------

  ArrayWithUnit<T, D> convert(const Unit& R) const
  { 
    ArrayWithUnit<T,D> res;
    res.value.reference(blitz::Array<T, D>(value * FullPhysics::conversion(units, R)));
    res.units = R;
    return res; 
  }

//-----------------------------------------------------------------------
/// We often need to handle conversion from wavenumber to/from
/// wavelength. This is either a normal conversion of the units before
/// and after match in the power of length (so cm^-1 to m^-1), or do
/// an inversion. Since we do this often enough, it is worth having a
/// function that handles this logic.
//-----------------------------------------------------------------------
    
  ArrayWithUnit<T,D> convert_wave(const Unit& R) const
  {
    ArrayWithUnit<T,D> res;
    res.units = R;
    if(units.is_commensurate(R))
      res.value.reference(blitz::Array<T, D>(FullPhysics::conversion(units, R) * value));
    else
      res.value.reference(blitz::Array<T, D>(FullPhysics::conversion(1 / units, R) / value));
    return res;
  }

  int rows() const {return value.rows();}
  int cols() const {return value.cols();}
  int depth() const {return value.depth();}

  void print(std::ostream& Os) const {
    Os << "ArrayWithUnit:\n"
       << "Value: " << value << "\n"
       << units << "\n";
  }

};
}
#endif
