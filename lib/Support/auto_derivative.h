#ifndef AUTO_DERIVATIVE_H
#define AUTO_DERIVATIVE_H
#include "fp_exception.h"
#include "printable.h"
#include <cmath>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>

// Because of the way blitz expands things, we need to declare these
// functions before including blitz.
namespace FullPhysics {
  template<class T> class AutoDerivative;
}
namespace std {			// Math functions are in std:: namespace.
FullPhysics::AutoDerivative<double> 
sqrt(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> 
log(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> 
log10(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> 
exp(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> 
sin(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> 
asin(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> 
cos(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> 
acos(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> 
tan(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> 
atan(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> 
pow(const FullPhysics::AutoDerivative<double>& x, const double y);
FullPhysics::AutoDerivative<double> 
pow(const double x, const FullPhysics::AutoDerivative<double>& y);
}
namespace blitz {
FullPhysics::AutoDerivative<double> 
pow2(const FullPhysics::AutoDerivative<double>& x);
}
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  Helper class that gives us a reference that we can assign a
  AutoDerivative to and write into the correct space in a ArrayAd.
*******************************************************************/
template<class T> class AutoDerivativeRef :
public Printable<AutoDerivativeRef<T> >
{
public:
  AutoDerivativeRef(T& V, const blitz::Array<T, 1>& G)
    : v(V), grad(G) 
  { }
  AutoDerivativeRef(T& V)
    : v(V)
  { }
  // This is const because although the values change, the reference doesn't
  const AutoDerivativeRef<T>& operator=(const AutoDerivative<T>& D) const
  {
    v = D.value();
    if(D.is_constant())
      grad = 0;
    else if(grad.rows() == D.gradient().rows())
      grad = D.gradient();
    else {
      Exception e;
      e << "Size of gradient in AutoDerivative doesn't match value you "
	<< "are trying to assign. Size of destination is " << grad.rows()
	<< " and size of source is " << D.gradient().rows();
      throw e;
    }
    return *this;
  }
  T& value_ref() const {return v;}
  blitz::Array<T, 1>& gradient_ref() const {return grad; }
  T value() const {return v;}
  const blitz::Array<T, 1>& gradient() const {return grad; }
  friend AutoDerivative<T> operator+(const AutoDerivativeRef<T>& X, 
				     const AutoDerivativeRef<T>& Y)
  { AutoDerivative<T> r(X);
    r+= Y;
    return r;
  }
  friend AutoDerivative<T> operator-(const AutoDerivativeRef<T>& X, 
				     const AutoDerivativeRef<T>& Y)
  { AutoDerivative<T> r(X);
    r-= Y;
    return r;
  }
  friend AutoDerivative<T> operator-(const AutoDerivativeRef<T>& X)
  { AutoDerivative<T> r(X);
    r*= -1;
    return r;
  }
  friend AutoDerivative<T> operator*(const AutoDerivativeRef<T>& X, 
				     const AutoDerivativeRef<T>& Y)
  { AutoDerivative<T> r(X);
    r*= Y;
    return r;
  }
  friend AutoDerivative<T> operator/(const AutoDerivativeRef<T>& X, 
				     const AutoDerivativeRef<T>& Y)
  { AutoDerivative<T> r(X);
    r/= Y;
    return r;
  }
  friend AutoDerivative<T> operator+(const AutoDerivativeRef<T>& X, 
				     const T& Y)
  { AutoDerivative<T> r(X);
    r+= Y;
    return r;
  }
  friend AutoDerivative<T> operator-(const AutoDerivativeRef<T>& X, 
				     const T& Y)
  { AutoDerivative<T> r(X);
    r-= Y;
    return r;
  }
  friend AutoDerivative<T> operator*(const AutoDerivativeRef<T>& X, 
				     const T& Y)
  { AutoDerivative<T> r(X);
    r*= Y;
    return r;
  }
  friend AutoDerivative<T> operator/(const AutoDerivativeRef<T>& X, 
				     const T& Y)
  { AutoDerivative<T> r(X);
    r/= Y;
    return r;
  }
  friend AutoDerivative<T> operator+(const T& X, 
				     const AutoDerivativeRef<T>& Y)
  { AutoDerivative<T> r(Y);
    r+= X;
    return r;
  }
  friend AutoDerivative<T> operator-(const T& X, 
				     const AutoDerivativeRef<T>& Y)
  { AutoDerivative<T> r(X);
    r-= Y;
    return r;
  }
  friend AutoDerivative<T> operator*(const T& X, 
				     const AutoDerivativeRef<T>& Y)
  { AutoDerivative<T> r(Y);
    r*= X;
    return r;
  }
  friend AutoDerivative<T> operator/(const T& X, 
				     const AutoDerivativeRef<T>& Y)
  { AutoDerivative<T> r(X);
    r/= Y;
    return r;
  }
  void print(std::ostream& Os) const
  { Os << "AutoDerivativeRef\n"
       << "  Value: " << value() << "\n"
       << "  Gradient:\n" << gradient() << "\n";
  }
private:
  T& v;
  // We can assign to grad because although the value changes the
  // reference doesn't
  mutable blitz::Array<T, 1> grad;
};

/****************************************************************//**
  There are a number of tools that can be used to do "Automatic
  Differentiation" (see for example
  http://www.autodiff.org/?module=Tools). 

  We examined several of the tools, and while these packages have a
  number of advantages (in particular, the ability to run the
  calculation either forward or backwards) for our particular needs a
  simpler forward only calculation was selected. This uses a number of
  type T, along with the first order gradient with respect to a set of
  independent variables.  We then overload the standard operations
  such as "+" and "*" to apply the chain rule, to propagate the
  derivatives forward.

  This is a newer field, so there doesn't seem to be standard
  terminology. In "Scientific and Engineering C++" by John Barton and
  Lee Nackman, this is call "Rall numbers" after a paper by L.B. Rall.
  On wikipedia at
  http://en.wikipedia.org/wiki/Automatic_differentiation this is
  called "Automatic differentiation using dual numbers". 

  As the existing automatic differentiation packages mature, we may
  want to revisit this choice and replace this class with a fuller
  library. 

  This class is not as efficient as hand coding derivative
  calculation, although it is much easier to use. If profiling shows a
  particular bottle neck, you might want to hand code derivatives for
  that specific location, which can then be placed into a
  AutoDerivative for use elsewhere (see the Absco class for an example
  of doing this).

  See also ArrayAd which work with this class for Arrays of 
  AutoDerivative. 
*******************************************************************/
template<class T> class AutoDerivative: 
public Printable<AutoDerivative<T> >,
boost::totally_ordered<class AutoDerivative<T> >,
boost::totally_ordered<class AutoDerivative<T>, T>
{
public:
  typedef T value_type;

//-----------------------------------------------------------------------
/// Default constructor, data is uninitialized.
//-----------------------------------------------------------------------

  AutoDerivative() {}

//-----------------------------------------------------------------------
/// Constructor that takes a value and a gradient. This makes a shallow
/// copy of the gradient, so we point to the same area in memory as G. If 
/// you want a deep copy, then pass in G.copy() rather than G.
//-----------------------------------------------------------------------
  
  AutoDerivative(const T& Val, const blitz::Array<T, 1>& G)
    : val(Val), grad(G) {}

//-----------------------------------------------------------------------
/// Convert AutoDerivativeRef to a AutoDerivative.
//-----------------------------------------------------------------------
 AutoDerivative(const AutoDerivativeRef<T>& V)
   : val(V.value()), grad(V.gradient().copy()) { }

//-----------------------------------------------------------------------
/// Constructor for a value of the i_th independent variable (0 based). We
/// create a gradient that is all 0, except for "1" in the i_th location.
//-----------------------------------------------------------------------

  AutoDerivative(const T& Val, int i_th, int nvars)
    : val(Val), grad(nvars)
  {
    range_check(i_th, 0, nvars);
    grad = 0;
    grad(i_th) = 1;
  }

//-----------------------------------------------------------------------
/// Constructor for a value that is a constant. 
//-----------------------------------------------------------------------

  AutoDerivative(const T& Val)
    : val(Val)
  {
  }

//-----------------------------------------------------------------------
/// Copy constructor. This does a deep copy.
//-----------------------------------------------------------------------

  AutoDerivative(const AutoDerivative<T>& D)
    : val(D.val), grad(D.grad.copy())
  {
  }

//-----------------------------------------------------------------------
/// Assignment operator. This makes a deep copy.
//-----------------------------------------------------------------------

  AutoDerivative<T>& operator=(const AutoDerivative<T>& D)
  {
    val = D.val;
    // Copy the data into existing gradient if it fits
    if(grad.rows() == D.grad.rows())
      grad = D.grad;
    else
      grad.reference(D.grad.copy());
    return *this;
  }

//-----------------------------------------------------------------------
/// Create deep copy.
//-----------------------------------------------------------------------

  AutoDerivative<T> copy() const
  { return AutoDerivative<T>(value(), gradient().copy()); }

//-----------------------------------------------------------------------
/// Convert to type T.
//-----------------------------------------------------------------------

  const T& value() const {return val;}
  T& value() {return val;}

//-----------------------------------------------------------------------
/// Gradient.
//-----------------------------------------------------------------------

  const blitz::Array<T, 1>& gradient() const { return grad; }
  blitz::Array<T, 1>& gradient() { return grad; }

//-----------------------------------------------------------------------
/// Number of variables in gradient.
//-----------------------------------------------------------------------

  int number_variable() const {return grad.rows();}

//-----------------------------------------------------------------------
/// Is this object a constant (with a gradient() all zeros)?
//-----------------------------------------------------------------------

  bool is_constant() const { return number_variable() == 0; }

//-----------------------------------------------------------------------
// Operators needed by boost to give all the arithmetic operators.
//-----------------------------------------------------------------------

  bool operator<(const AutoDerivative<T>& V) const {return val < V.val;}
  bool operator<(const T& V) const {return val < V;}
  bool operator>(const AutoDerivative<T>& V) const {return val > V.val;}
  bool operator>(const T& V) const {return val > V;}
  bool operator==(const AutoDerivative<T>& V) const {return val == V.val;}
  bool operator==(const T& V) const {return val == V;}
  AutoDerivative<T>& operator+=(const AutoDerivative<T>& V)
  { 
    val += V.val; 
    // Handling for V or self actually being a constant, which means
    // we have a zero size gradient.
    if(!V.is_constant()) {
      if(!is_constant())
	grad += V.grad; 
      else
	grad.reference(V.grad.copy());
    }
    return *this;
  }
  AutoDerivative<T>& operator+=(const T& V)
  { val += V; return *this;}
  AutoDerivative<T>& operator-=(const AutoDerivative<T>& V)
  { val -= V.val; 
    // Handling for V or self actually being a constant, which means
    // we have a zero size gradient.
    if(!V.is_constant()) {
      if(!is_constant())
	grad -= V.grad; 
      else {
	grad.resize(V.grad.shape());
	grad = -V.grad;
      }
    }
    return *this;
  }
  AutoDerivative<T>& operator-=(const T& V)
  { val -= V; return *this;}
  AutoDerivative<T>& operator*=(const AutoDerivative<T>& V)
  { 
    // Handling for V or self actually being a constant, which means
    // we have a zero size gradient.
    if(!V.is_constant()) {
      if(!is_constant())
	grad = V.val * grad + val * V.grad; 
      else {
	grad.reference(V.grad.copy());
	grad *= val;
      }
    } else
      grad *= V.val; 
    val *= V.val; 
    return *this;
  }
  AutoDerivative<T>& operator*=(const T& V)
  { val *= V; grad *= V; return *this;}
  AutoDerivative<T>& operator/=(const AutoDerivative<T>& V)
  { 
    // Handling for V or self actually being a constant, which means
    // we have a zero size gradient.
    if(!V.is_constant()) {
      if(!is_constant())
	grad = 1 / V.val * grad - val / (V.val * V.val) * V.grad; 
      else {
	grad.resize(V.grad.shape());
	grad = - val / (V.val * V.val) * V.grad; 
      }
    } else
      grad /= V.val; 
    val /= V.val; 
    return *this;
  }
  AutoDerivative<T>& operator/=(const T& V)
  { val /= V; grad /= V; return *this;}
  void print(std::ostream& Os) const
  { Os << "AutoDerivative\n"
       << "  Value: " << val << "\n"
       << "  Gradient:\n" << grad << "\n";
  }
  friend AutoDerivative<T> operator+(const AutoDerivative<T>& X, 
				     const AutoDerivative<T>& Y)
  { AutoDerivative<T> r(X);
    r+= Y;
    return r;
  }
  friend AutoDerivative<T> operator-(const AutoDerivative<T>& X, 
				     const AutoDerivative<T>& Y)
  { AutoDerivative<T> r(X);
    r -= Y;
    return r;
  }
  friend AutoDerivative<T> operator-(const AutoDerivative<T>& X)
  { AutoDerivative<T> r(X);
    r *= -1;
    return r;
  }
  friend AutoDerivative<T> operator*(const AutoDerivative<T>& X, 
				     const AutoDerivative<T>& Y)
  { AutoDerivative<T> r(X);
    r *= Y;
    return r;
  }
  friend AutoDerivative<T> operator/(const AutoDerivative<T>& X, 
				     const AutoDerivative<T>& Y)
  { AutoDerivative<T> r(X);
    r /= Y;
    return r;
  }
  friend AutoDerivative<T> operator+(const AutoDerivative<T>& X, 
				     const T& Y)
  { AutoDerivative<T> r(X);
    r+= Y;
    return r;
  }
  friend AutoDerivative<T> operator-(const AutoDerivative<T>& X, 
				     const T& Y)
  { AutoDerivative<T> r(X);
    r -= Y;
    return r;
  }
  friend AutoDerivative<T> operator*(const AutoDerivative<T>& X, 
				     const T& Y)
  { AutoDerivative<T> r(X);
    r *= Y;
    return r;
  }
  friend AutoDerivative<T> operator/(const AutoDerivative<T>& X, 
				     const T& Y)
  { AutoDerivative<T> r(X);
    r /= Y;
    return r;
  }
  friend AutoDerivative<T> operator+(const T& X, 
				     const AutoDerivative<T>& Y)
  { AutoDerivative<T> r(Y);
    r+= X;
    return r;
  }
  friend AutoDerivative<T> operator-(const T& X, 
				     const AutoDerivative<T>& Y)
  { AutoDerivative<T> r(X);
    r -= Y;
    return r;
  }
  friend AutoDerivative<T> operator*(const T& X, 
				     const AutoDerivative<T>& Y)
  { AutoDerivative<T> r(Y);
    r *= X;
    return r;
  }
  friend AutoDerivative<T> operator/(const T& X, 
				     const AutoDerivative<T>& Y)
  { AutoDerivative<T> r(X);
    r /= Y;
    return r;
  }
private:
  T val;
  blitz::Array<T, 1> grad;
  bool keep_grad;
};

//-----------------------------------------------------------------------
/// Utility function to extract the Jacobian as a separate matrix from
/// an array of AutoDerivative.
//-----------------------------------------------------------------------

template<class T> inline blitz::Array<T, 2> jacobian
(const blitz::Array<AutoDerivative<T>, 1>& Ad)
{
  int nvar = 0;
  BOOST_FOREACH(const AutoDerivative<T>& v, Ad)
    if(!v.is_constant()) {
      nvar = v.number_variable();
      break;
    }
  blitz::Array<T, 2> res(Ad.rows(), nvar);
  if(nvar > 0) {
    for(int i = 0; i < Ad.rows(); ++i)
      if(!Ad(i).is_constant())
	res(i, blitz::Range::all()) = Ad(i).gradient();
      else
	res(i, blitz::Range::all()) = 0;
  }
  return res;
}

template<class T> inline blitz::Array<T, 3> jacobian
(const blitz::Array<AutoDerivative<T>, 2>& Ad)
{
  int nvar = 0;
  BOOST_FOREACH(const AutoDerivative<T>& v, Ad)
    if(!v.is_constant()) {
      nvar = v.number_variable();
      break;
    }
  blitz::Array<T, 3> res(Ad.rows(), Ad.cols(), nvar);
  if(nvar > 0) {
    for(int i = 0; i < Ad.rows(); ++i)
      for(int j = 0; j < Ad.cols(); ++j)
	if(!Ad(i, j).is_constant())
	  res(i, j, blitz::Range::all()) = Ad(i, j).gradient();
	else
	  res(i, j, blitz::Range::all()) = 0;
  }
  return res;
}

template<class T> inline blitz::Array<T, 4> jacobian
(const blitz::Array<AutoDerivative<T>, 3>& Ad)
{
  int nvar = 0;
  BOOST_FOREACH(const AutoDerivative<T>& v, Ad)
    if(!v.is_constant()) {
      nvar = v.number_variable();
      break;
    }
  blitz::Array<T, 4> res(Ad.rows(), Ad.cols(), Ad.depth(), nvar);
  if(nvar > 0) {
    for(int i = 0; i < Ad.rows(); ++i)
      for(int j = 0; j < Ad.cols(); ++j)
	for(int k = 0; k < Ad.depth(); ++k)
	  if(!Ad(i, j, k).is_constant())
	    res(i, j, k, blitz::Range::all()) = Ad(i, j, k).gradient();
	  else
	    res(i, j, k, blitz::Range::all()) = 0;
  }
  return res;
}

template<class T> inline blitz::Array<T, 5> jacobian
(const blitz::Array<AutoDerivative<T>, 4>& Ad)
{
  int nvar = 0;
  BOOST_FOREACH(const AutoDerivative<T>& v, Ad)
    if(!v.is_constant()) {
      nvar = v.number_variable();
      break;
    }
  blitz::Array<T, 5> res(Ad.rows(), Ad.cols(), Ad.depth(), 
			 Ad.extent(blitz::fourthDim), nvar);
  if(nvar > 0) {
    for(int i = 0; i < Ad.rows(); ++i)
      for(int j = 0; j < Ad.cols(); ++j)
	for(int k = 0; k < Ad.depth(); ++k)
	  for(int m = 0; m < Ad.extent(blitz::fourthDim); ++m)
	    if(!Ad(i, j, k, m).is_constant())
	      res(i, j, k, m, blitz::Range::all()) = Ad(i, j, k, m).gradient();
	    else
	      res(i, j, k, m, blitz::Range::all()) = 0;
  }
  return res;
}

//-----------------------------------------------------------------------
/// Utility routine to take value and Jacobian separately and create
/// a blitz::Array<AutoDerivative<double>, 1>
//-----------------------------------------------------------------------

inline blitz::Array<AutoDerivative<double>, 1> auto_derivative
(const blitz::Array<double, 1>& Val, const blitz::Array<double, 2>& Jac)
{
  blitz::Array<AutoDerivative<double>, 1> res(Val.shape());
  for(int i = 0; i < Val.rows(); ++i) {
    res(i).value() = Val(i);
    if(Jac.cols() > 0)
      res(i).gradient().reference(Jac(i, blitz::Range::all()));
  }
  return res;
}

inline blitz::Array<AutoDerivative<double>, 1> auto_derivative
(const blitz::Array<double, 1>& Val, int i_th, int nvars)
{
  blitz::Array<AutoDerivative<double>, 1> res(Val.shape());
  for(int i = 0; i < Val.rows(); ++i)
    res(i) = AutoDerivative<double>(Val(i), i_th, nvars);
  return res;
}

inline blitz::Array<AutoDerivative<double>, 2> auto_derivative
(const blitz::Array<double, 2>& Val, const blitz::Array<double, 3>& Jac)
{
  blitz::Array<AutoDerivative<double>, 2> res(Val.shape());
  for(int i = 0; i < Val.rows(); ++i)
    for(int j = 0; j < Val.cols(); ++j) {
      res(i, j).value() = Val(i, j);
      if(Jac.depth() > 0)
	res(i, j).gradient().reference(Jac(i, j, blitz::Range::all()));
    }
  return res;
}

inline blitz::Array<AutoDerivative<double>, 3> auto_derivative
(const blitz::Array<double, 3>& Val, const blitz::Array<double, 4>& Jac)
{
  blitz::Array<AutoDerivative<double>, 3> res(Val.shape());
  for(int i = 0; i < Val.rows(); ++i)
    for(int j = 0; j < Val.cols(); ++j)
      for(int k = 0; k < Val.extent(blitz::thirdDim); ++k) {
	res(i, j, k).value() = Val(i, j, k);
	if(Jac.extent(blitz::fourthDim) > 0)
	  res(i, j, k).gradient().reference(Jac(i, j, k, blitz::Range::all()));
      }
  return res;
}
}

//-----------------------------------------------------------------------
/// Math functions.
//-----------------------------------------------------------------------
namespace std {			// Math functions are in std:: namespace.

inline FullPhysics::AutoDerivative<double> 
sqrt(const FullPhysics::AutoDerivative<double>& x)
{
  if(x.is_constant())
    return FullPhysics::AutoDerivative<double>(::sqrt(x.value()));
  else
    return FullPhysics::AutoDerivative<double>(::sqrt(x.value()), 
     blitz::Array<double, 1>(x.gradient() / (2*::sqrt(x.value())) ));
}

inline FullPhysics::AutoDerivative<double> 
log(const FullPhysics::AutoDerivative<double>& x)
{
  if(x.is_constant())
    return FullPhysics::AutoDerivative<double>(::log(x.value()));
  else
    return FullPhysics::AutoDerivative<double>(::log(x.value()), 
       blitz::Array<double, 1>(x.gradient() / x.value()));
}

inline FullPhysics::AutoDerivative<double> 
log10(const FullPhysics::AutoDerivative<double>& x)
{
  if(x.is_constant())
    return FullPhysics::AutoDerivative<double>(::log10(x.value()));
  else
    return FullPhysics::AutoDerivative<double>(::log10(x.value()), 
       blitz::Array<double, 1>(M_LOG10E * x.gradient() / x.value()));
}

inline FullPhysics::AutoDerivative<double> 
exp(const FullPhysics::AutoDerivative<double>& x)
{
  if(x.is_constant())
    return FullPhysics::AutoDerivative<double>(::exp(x.value()));
  else 
    return FullPhysics::AutoDerivative<double>(::exp(x.value()), 
     blitz::Array<double, 1>(x.gradient() * ::exp(x.value())));
}

inline FullPhysics::AutoDerivative<double> 
sin(const FullPhysics::AutoDerivative<double>& x)
{
  if(x.is_constant())
    return FullPhysics::AutoDerivative<double>(::sin(x.value()));
  else {
    return FullPhysics::AutoDerivative<double>(::sin(x.value()),
     blitz::Array<double, 1>(x.gradient() * ::cos(x.value()) ));
  }
}

inline FullPhysics::AutoDerivative<double> 
asin(const FullPhysics::AutoDerivative<double>& x)
{
  if(x.value() < -1.0 || x.value() > 1.0) {
    FullPhysics::Exception err_msg;
    err_msg << "Value not within domain of arcsine: -1 <= x <= 1 where x = " << x.value();
    throw err_msg;
  }
    
  if(x.is_constant())
    return FullPhysics::AutoDerivative<double>(::asin(x.value()));
  else {
    return FullPhysics::AutoDerivative<double>(::asin(x.value()),
     blitz::Array<double, 1>(x.gradient() / sqrt(1-x.value()*x.value()) ));
  }
}

inline FullPhysics::AutoDerivative<double> 
cos(const FullPhysics::AutoDerivative<double>& x)
{
  if(x.is_constant())
    return FullPhysics::AutoDerivative<double>(::cos(x.value()));
  else {
    return FullPhysics::AutoDerivative<double>(::cos(x.value()),
     blitz::Array<double, 1>(x.gradient() * (- ::sin(x.value())) ));
  }
}

inline FullPhysics::AutoDerivative<double> 
acos(const FullPhysics::AutoDerivative<double>& x)
{
  if(x.value() < -1.0 || x.value() > 1.0) {
    FullPhysics::Exception err_msg;
    err_msg << "Value not within domain of arccosine: -1 <= x <= 1 where x = " << x.value();
    throw err_msg;
  }
    
  if(x.is_constant())
    return FullPhysics::AutoDerivative<double>(::acos(x.value()));
  else {
    return FullPhysics::AutoDerivative<double>(::acos(x.value()),
     blitz::Array<double, 1>(-1.0 * x.gradient() / sqrt(1-x.value()*x.value()) ));
  }
}

inline FullPhysics::AutoDerivative<double> 
tan(const FullPhysics::AutoDerivative<double>& x)
{
  if(x.is_constant())
    return FullPhysics::AutoDerivative<double>(::tan(x.value()));
  else {
    return FullPhysics::AutoDerivative<double>(::tan(x.value()),
     blitz::Array<double, 1>(x.gradient() / (::cos(x.value()) * ::cos(x.value())) ));
  }
}

inline FullPhysics::AutoDerivative<double> 
atan(const FullPhysics::AutoDerivative<double>& x)
{   
  if(x.is_constant())
    return FullPhysics::AutoDerivative<double>(::atan(x.value()));
  else {
    return FullPhysics::AutoDerivative<double>(::atan(x.value()),
     blitz::Array<double, 1>(x.gradient() / (x.value()*x.value() + 1) ));
  }
}

inline FullPhysics::AutoDerivative<double> 
pow(const FullPhysics::AutoDerivative<double>& base, const double exponent)
{   
  if(base.is_constant())
    return FullPhysics::AutoDerivative<double>(::pow(base.value(), exponent));
  else {
    return FullPhysics::AutoDerivative<double>(::pow(base.value(), exponent),
      blitz::Array<double, 1>( base.gradient() * exponent * pow(base.value(), exponent - 1) ));
  }
}

inline FullPhysics::AutoDerivative<double> 
pow(const double base, const FullPhysics::AutoDerivative<double>& exponent)
{   
  if(exponent.is_constant())
    return FullPhysics::AutoDerivative<double>(::pow(base, exponent.value()));
  else {
    return FullPhysics::AutoDerivative<double>(::pow(base, exponent.value()),
       blitz::Array<double, 1>( exponent.gradient() * ::log(base) * ::pow(base, exponent.value()) ));
  }
}

}


//-----------------------------------------------------------------------
/// Apply value function to a blitz array.
//-----------------------------------------------------------------------

namespace blitz {
  inline double value(const FullPhysics::AutoDerivative<double>& Ad )
  { return Ad.value(); }
  BZ_DECLARE_FUNCTION_RET(value, double);
  inline FullPhysics::AutoDerivative<double> auto_derivative(double Ad) 
  { return FullPhysics::AutoDerivative<double>(Ad); }
  BZ_DECLARE_FUNCTION_RET(auto_derivative, FullPhysics::AutoDerivative<double>);
  BZ_DECLARE_ARRAY_ET_SCALAR_FUNCS(FullPhysics::AutoDerivative<double>)
  BZ_DECLARE_ARRAY_ET_SCALAR_OPS(FullPhysics::AutoDerivative<double>)

inline FullPhysics::AutoDerivative<double> 
pow2(const FullPhysics::AutoDerivative<double>& base)
{   
  if(base.is_constant())
    return FullPhysics::AutoDerivative<double>(base.value() * base.value());
  else {
    return FullPhysics::AutoDerivative<double>(base.value() * base.value(),
       blitz::Array<double, 1>( base.gradient() * 2 * base.value()));
  }
}
}
#endif
