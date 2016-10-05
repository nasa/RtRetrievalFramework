#ifndef ARRAY_AD_H
#define ARRAY_AD_H
#include "auto_derivative.h"
#include "fp_exception.h"

// Turn on trace messages, useful for debugging problems with this class.
// #define ARRAY_AD_DIAGNOSTIC
#ifdef ARRAY_AD_DIAGNOSTIC
#define ARRAY_AD_DIAGNOSTIC_MSG std::cerr << "In " << __LINE__ << "\n"
#else
#define ARRAY_AD_DIAGNOSTIC_MSG 
#endif

namespace FullPhysics {
/****************************************************************//**
  The AutoDerivative<T> works well, and it works with blitz if you 
  create a blitz::Array<AutoDerivative<T>, D>, however the performance
  tends to be a bit poor for larger Arrays because of all the 
  temporaries that get created.

  This class keeps all the values together in one blitz::Array, and all 
  the gradients in one Jacobian blitz::Array. This is less flexible than
  a blitz::Array<AutoDerivative<T>, D>, but the peformance is much better.
  We can convert to and from AutoDerivative<T> as needed, so this gives us
  the best of both worlds - ArrayAd when we don't need the extra 
  functionality, convert to AutoDerivative<T> when needed.

  A note on the implementation, blitz doesn't handle Range::all() well with 
  a zero size array. So we always ensure that the Jacobian has at least 
  one column, even for constant values with the nvars = 0. We track if
  this is a constant separately.

  We can get odd floating point errors (particularly on Ubuntu), 
  if we copy garbage values. So even though we don't use the
  Jacobian, if the data is constant, we initialize it to 0.
  
*******************************************************************/

template<class T, int D> class ArrayAd 
{
public:
  ArrayAd(const blitz::Array<AutoDerivative<T>, D>& V)
    : val(V.shape()), jac(FullPhysics::jacobian(V)), 
      is_const(false) 
  { 
    ARRAY_AD_DIAGNOSTIC_MSG;
    if(number_variable() ==0)
      resize(V.shape(), 0);
    val = blitz::value(V);
  }
  ArrayAd(const ArrayAd<T, D>& V)
    : val(V.val), jac(V.jac), is_const(V.is_const) 
  { 
    ARRAY_AD_DIAGNOSTIC_MSG;
  }
  ArrayAd& operator=(const ArrayAd<T, D>& V)
  { 
    ARRAY_AD_DIAGNOSTIC_MSG;
    if(val.rows() ==0) { // Special handling for new empty objects
      val.reference(V.val.copy());
      jac.reference(V.jac.copy());
      is_const = V.is_const;
    } else {
      val = V.val; 
      if(V.is_const) {
	if(!is_const)
	  jac = 0;
      } else {
      // This is such a common error, add a check here.
	if(number_variable() != V.number_variable()) {
	  Exception e;
	  e << "Trying to assign Jacobian with " << V.number_variable() 
	    << " variables to size " << number_variable() << "\n";
	  throw e;
	}
	jac = V.jac;
      }
    }
    return *this;
  }
  ArrayAd(const blitz::Array<T, D>& Val, const blitz::Array<T, D+1>& Jac, 
	  bool Is_const = false, bool Force_copy = false)
    : val(Force_copy ? Val.copy() : Val), 
      jac(Force_copy ? Jac.copy() : Jac), is_const(Is_const) 
  {
    ARRAY_AD_DIAGNOSTIC_MSG;
    if(jac.extent(D) == 0) {
      is_const = true;
      blitz::TinyVector<int, D + 1> s2;
      for(int i = 0; i < D; ++i)
	s2(i) = val.shape()(i);
      s2(D) = 1;
      jac.resize(s2);
      jac = 0;
    }
  }
  ArrayAd(const blitz::Array<T, D>& Val, bool Force_copy = false)
    : val(Force_copy ? Val.copy() : Val), is_const(true) 
  {
    ARRAY_AD_DIAGNOSTIC_MSG;
    blitz::TinyVector<int, D + 1> s2;
    for(int i = 0; i < D; ++i)
      s2(i) = val.shape()(i);
    s2(D) = 1;
    jac.resize(s2);
    jac = 0;
  }
  ArrayAd() :val(0), jac(0,1), is_const(false) {}
  ArrayAd(int n1, int nvar)
    : val(n1), jac(n1, std::max(nvar, 1)), is_const(nvar == 0) 
  {
    if(is_const)
      jac = 0;
  }
  ArrayAd(int n1, int n2, int nvar)
    : val(n1, n2), jac(n1, n2, std::max(nvar, 1)), is_const(nvar == 0) 
  {
    if(is_const)
      jac = 0;
  }
  ArrayAd(int n1, int n2, int n3, int nvar)
    : val(n1, n2, n3), jac(n1, n2, n3, std::max(nvar, 1)), 
      is_const(nvar == 0) 
  {
    if(is_const)
      jac = 0;
  }
  ArrayAd(int n1, int n2, int n3, int n4, int nvar)
    : val(n1, n2, n3, n4), jac(n1, n2, n3, n4, std::max(nvar, 1)), 
      is_const(nvar == 0) 
  {
    if(is_const)
      jac = 0;
  }
  ArrayAd(int n1, int n2, int n3, int n4, int n5, int nvar)
    : val(n1, n2, n3, n4, n5), jac(n1, n2, n3, n4, n5, std::max(nvar, 1)), 
      is_const(nvar == 0) 
  {
    if(is_const)
      jac = 0;
  }
  ArrayAd(const blitz::TinyVector<int, D>& Shape, int nvar)
    : val(Shape), is_const(nvar == 0)
  {
    blitz::TinyVector<int, D + 1> s2;
    for(int i = 0; i < D; ++i)
      s2(i) = Shape(i);
    s2(D) = std::max(nvar, 1);
    jac.resize(s2);
    if(is_const)
      jac = 0;
  }
  void resize_number_variable(int nvar)
  {
    if(nvar == number_variable())
      return;
    blitz::TinyVector<int, D + 1> s2;
    for(int i = 0; i < D; ++i)
      s2(i) = val.shape()(i);
    s2(D) = std::max(nvar, 1);
    jac.resize(s2);
    jac = 0;
    is_const = (nvar ==0);
  }
  void resize(const blitz::TinyVector<int, D>& Shape, int nvar)
  {
    val.resize(Shape);
    is_const = (nvar ==0);
    blitz::TinyVector<int, D + 1> s2;
    for(int i = 0; i < D; ++i)
      s2(i) = Shape(i);
    s2(D) = std::max(nvar, 1);
    jac.resize(s2);
    if(is_const)
      jac = 0;
  }
  void resize(int n1, int nvar)
  { 
    val.resize(n1); jac.resize(n1, std::max(nvar, 1)); 
    is_const = (nvar == 0); 
    if(is_const)
      jac = 0;
  }
  void resize(int n1, int n2, int nvar)
  { 
    val.resize(n1, n2); jac.resize(n1, n2, std::max(nvar, 1)); 
    is_const = (nvar == 0); 
    if(is_const)
      jac = 0;
  }
  void resize(int n1, int n2, int n3, int nvar)
  { 
    val.resize(n1, n2, n3); jac.resize(n1, n2, n3, std::max(nvar, 1)); 
    is_const = (nvar == 0); 
    if(is_const)
      jac = 0;
  }
  void resize(int n1, int n2, int n3, int n4, int nvar)
  { 
    val.resize(n1, n2, n3, n4); 
    jac.resize(n1, n2, n3, n4, std::max(nvar, 1)); 
    is_const = (nvar == 0); 
    if(is_const)
      jac = 0;
  }
  void resize(int n1, int n2, int n3, int n4, int n5, int nvar)
  { 
    val.resize(n1, n2, n3, n4, n5); 
    jac.resize(n1, n2, n3, n4, n5, std::max(nvar, 1)); 
    is_const = (nvar == 0); 
    if(is_const)
      jac = 0;
  }
  AutoDerivativeRef<T> operator()(int i1)
  { 
    if(is_const)
      return AutoDerivativeRef<T>(val(i1));
    else
      return AutoDerivativeRef<T>(val(i1), 
				  jac(i1, blitz::Range::all())); 
  }
  AutoDerivativeRef<T> operator()(int i1, int i2)
  { 
    if(is_const)
      return AutoDerivativeRef<T>(val(i1, i2));
    else
      return AutoDerivativeRef<T>(val(i1, i2),
				  jac(i1, i2, blitz::Range::all())); 
  }
  AutoDerivativeRef<T> operator()(int i1, int i2, int i3)
  { 
    if(is_const)
      return AutoDerivativeRef<T>(val(i1, i2, i3));
    else
      return AutoDerivativeRef<T>(val(i1, i2, i3), 
				  jac(i1, i2, i3, blitz::Range::all())); 
  }
  AutoDerivativeRef<T> operator()(int i1, int i2, int i3, int i4)
  { 
    if(is_const)
      return AutoDerivativeRef<T>(val(i1, i2, i3, i4));
    else
      return AutoDerivativeRef<T>(val(i1, i2, i3, i4), 
				  jac(i1, i2, i3, i4, blitz::Range::all())); 
  }
  AutoDerivativeRef<T> operator()(int i1, int i2, int i3, int i4, int i5)
  { 
    if(is_const)
      return AutoDerivativeRef<T>(val(i1, i2, i3, i4, i5));
    else
      return AutoDerivativeRef<T>(val(i1, i2, i3, i4, i5), 
				  jac(i1, i2, i3, i4, i5, blitz::Range::all())); 
  }
  AutoDerivative<T> operator()(int i1) const
  { 
    if(is_const)
      return AutoDerivative<T>(val(i1));
    else
      return AutoDerivative<T>(val(i1), 
			       jac(i1, blitz::Range::all())); 
  }
  AutoDerivative<T> operator()(int i1, int i2) const
  { 
    if(is_const)
      return AutoDerivative<T>(val(i1, i2));
    else
      return AutoDerivative<T>(val(i1, i2),
			       jac(i1, i2, blitz::Range::all())); 
  }
  AutoDerivative<T> operator()(int i1, int i2, int i3) const
  { 
    if(is_const)
      return AutoDerivative<T>(val(i1, i2, i3));
    else
      return AutoDerivative<T>(val(i1, i2, i3), 
			       jac(i1, i2, i3, blitz::Range::all())); 
  }
  AutoDerivative<T> operator()(int i1, int i2, int i3, int i4) const
  { 
    if(is_const)
      return AutoDerivative<T>(val(i1, i2, i3, i4));
    else
      return AutoDerivative<T>(val(i1, i2, i3, i4), 
			       jac(i1, i2, i3, i4, blitz::Range::all())); 
  }
  AutoDerivative<T> operator()(int i1, int i2, int i3, int i4, int i5) const
  { 
    if(is_const)
      return AutoDerivative<T>(val(i1, i2, i3, i4, i5));
    else
      return AutoDerivative<T>(val(i1, i2, i3, i4, i5), 
			       jac(i1, i2, i3, i4, i5, blitz::Range::all())); 
  }
  const blitz::Array<T, D>& value() const {return val;}
  const blitz::Array<T, D+1> jacobian() const 
  {
    if(is_const)
      return blitz::Array<T, D+1>();
    return jac;
  }
  blitz::Array<T, D>& value() {return val;}
  blitz::Array<T, D+1>& jacobian() 
  { return jac;}
  blitz::Array<AutoDerivative<T>, D> to_array() const
  { if(!is_const)
      return auto_derivative(val, jac); 
    blitz::Array<AutoDerivative<T>, D> res(val.shape());
    res = auto_derivative(val);
    return res;
  }
  ArrayAd& operator=(const blitz::Array<T, D>& V)
  {  ARRAY_AD_DIAGNOSTIC_MSG; val = V; jac = 0; return *this; }
  ArrayAd& operator=(const T& V)
  { ARRAY_AD_DIAGNOSTIC_MSG; val = V; jac = 0; return *this; }
  ArrayAd& operator=(const AutoDerivative<T>& V)
  { ARRAY_AD_DIAGNOSTIC_MSG; 
    val = V.value(); 
    blitz::IndexPlaceholder<D> ig;
    if(!is_constant())
      jac = V.gradient()(ig); 
    return *this; 
  }
  ArrayAd<T, 1> operator()(const blitz::Range& r1) const
  {
    return ArrayAd<T, 1>(val(r1), jac(r1, blitz::Range::all()), is_const);
  }
  template<typename T1, typename T2>
  ArrayAd<T,blitz::SliceInfo<T,T1,T2>::rank> operator()(T1 r1, T2 r2) const
  {
    return ArrayAd<T, blitz::SliceInfo<T,T1,T2>::rank>
      (val(r1, r2), jac(r1, r2, blitz::Range::all()), is_const);
  }
  template<typename T1, typename T2, typename T3>
  ArrayAd<T,blitz::SliceInfo<T,T1,T2,T3>::rank> 
  operator()(T1 r1, T2 r2, T3 r3) const
  {
    return ArrayAd<T, blitz::SliceInfo<T,T1,T2,T3>::rank>
      (val(r1, r2, r3), jac(r1, r2, r3, blitz::Range::all()), is_const);
  }
  template<typename T1, typename T2, typename T3, typename T4>
  ArrayAd<T,blitz::SliceInfo<T,T1,T2,T3,T4>::rank> 
  operator()(T1 r1, T2 r2, T3 r3, T4 r4) const
  {
    return ArrayAd<T, blitz::SliceInfo<T,T1,T2,T3,T4>::rank>
      (val(r1, r2, r3, r4), jac(r1, r2, r3, r4, blitz::Range::all()),
       is_const);
  }
  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  ArrayAd<T,blitz::SliceInfo<T,T1,T2,T3,T4,T5>::rank> 
  operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5) const
  {
    return ArrayAd<T, blitz::SliceInfo<T,T1,T2,T3,T4,T5>::rank>
      (val(r1, r2, r3, r4, r5), jac(r1, r2, r3, r4, r5, blitz::Range::all()),
       is_const);
  }
  int rows() const {return val.rows();}
  int cols() const {return val.cols();}
  int depth() const {return val.depth();}
  bool is_constant() const {return is_const;}
  void reference(const ArrayAd<T, D>& V)
  { val.reference(V.val); jac.reference(V.jac); is_const = V.is_const;}
  ArrayAd<T, D> copy() const
  { return ArrayAd<T,D>(val.copy(), jac.copy(), is_const); }
  int number_variable() const 
  {
    if(is_const)
      return 0;
    else
      return jac.extent(D);
  }
  const blitz::TinyVector<int, D>& shape() const { return val.shape(); }
  friend std::ostream& operator<<(std::ostream& os, 
				  const ArrayAd<T, D>& V)
  { os << V.val << "\n" << V.jac << "\n"; return os;}
  friend std::istream& operator>>(std::istream& is, ArrayAd<T, D>& V)
  { is >> V.val >> V.jac; return is; }
private:
  blitz::Array<T, D> val;
  blitz::Array<T, D + 1> jac;
  bool is_const;
};
}
#endif

