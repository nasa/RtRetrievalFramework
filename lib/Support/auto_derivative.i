// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "std_vector.i"
%include "common.i"

%{
#include "auto_derivative.h"
namespace FullPhysics {
}

%}
%base_import(generic_object)
namespace FullPhysics {
template<class T> class AutoDerivative;
template<class T> class AutoDerivativeRef : public GenericObject
{
public:
  AutoDerivativeRef(T& V, blitz::Array<T,1>& FORCE_COPY);
  %python_attribute(value, T)
  %python_attribute(gradient, blitz::Array<T, 1>)
  std::string print_to_string() const;
};

template<class T> class AutoDerivative : public GenericObject
{
public:
  typedef T value_type;
  AutoDerivative();
  AutoDerivative(const T& Val, const blitz::Array<T, 1>& FORCE_COPY);
  AutoDerivative(const T& Val, int i_th, int nvars);
  AutoDerivative(const T& Val);
  AutoDerivative(const AutoDerivative<T>& D);
  AutoDerivative(const AutoDerivativeRef<T>& V);
  %python_attribute(number_variable, int)
  %python_attribute(is_constant, bool)
  bool operator<(const AutoDerivative<T>& V) const;
  bool operator<(const T& V) const;
  bool operator>(const AutoDerivative<T>& V) const;
  bool operator>(const T& V) const;
  bool operator==(const AutoDerivative<T>& V) const;
  bool operator==(const T& V) const;
  AutoDerivative<T> operator+=(const AutoDerivative<T>& V);
  AutoDerivative<T> operator+=(const T& V);
  AutoDerivative<T> operator-=(const AutoDerivative<T>& V);
  AutoDerivative<T> operator-=(const T& V);
  AutoDerivative<T> operator*=(const AutoDerivative<T>& V);
  AutoDerivative<T> operator*=(const T& V);
  AutoDerivative<T> operator/=(const AutoDerivative<T>& V);
  AutoDerivative<T> operator/=(const T& V);
  std::string print_to_string() const;
  %extend {
    T _value() const { return $self->value();}
    void _value_set(T V) { $self->value() = V;}
    blitz::Array<T, 1> _gradient() const { return $self->gradient();}
    void _gradient_set(blitz::Array<T, 1>& V) { $self->gradient().reference(V.copy());}
    AutoDerivative<T> __add__(const AutoDerivative<T>& Y) 
    { return *$self + Y; }
    AutoDerivative<T> __add__(const T& Y) 
    { return *$self + Y; }
    AutoDerivative<T> __radd__(const T& X) 
    { return X + *$self;}
    AutoDerivative<T> __sub__(const AutoDerivative<T>& Y) 
    { return *$self - Y; }
    AutoDerivative<T> __sub__(const T& Y) 
    { return *$self - Y; }
    AutoDerivative<T> __rsub__(const T& X) 
    { return X - *$self;}
    AutoDerivative<T> __mul__(const AutoDerivative<T>& Y) 
    { return *$self * Y; }
    AutoDerivative<T> __mul__(const T& Y) 
    { return *$self * Y; }
    AutoDerivative<T> __rmul__(const T& X) 
    { return X * *$self;}
    // Python 2 division operator name
    AutoDerivative<T> __div__(const AutoDerivative<T>& Y) 
    { return *$self / Y; }
    AutoDerivative<T> __div__(const T& Y) 
    { return *$self / Y; }
    // Python 3 division operator name
    AutoDerivative<T> __truediv__(const AutoDerivative<T>& Y) 
    { return *$self / Y; }
    AutoDerivative<T> __truediv__(const T& Y) 
    { return *$self / Y; }
     AutoDerivative<T> __rdiv__(const T& X) 
    { return X / *$self;}
    AutoDerivative<T> __pow__(const T& X) 
    { return std::pow(*$self,  X);}
    AutoDerivative<T> __rpow__(const T& X) 
    { return std::pow(X, *$self);}
  }
  %pythoncode {
@property
def value(self):
  return self._value()

@value.setter
def value(self, val):
  self._value_set(val)

@property
def gradient(self):
  return self._gradient()

@gradient.setter
def gradient(self,val):
    self._gradient_set(val)
  }

};
}

%fp_shared_ptr(FullPhysics::AutoDerivative<double>)
%fp_shared_ptr(FullPhysics::AutoDerivativeRef<double>)
%template(AutoDerivativeDouble) FullPhysics::AutoDerivative<double>;
%template(AutoDerivativeRefDouble) FullPhysics::AutoDerivativeRef<double>;
%template(ArrayAutoDerivative_double_1) blitz::Array<FullPhysics::AutoDerivative<double>, 1>;

namespace std {
FullPhysics::AutoDerivative<double> sqrt(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> log(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> log10(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> exp(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> sin(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> asin(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> cos(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> acos(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> tan(const FullPhysics::AutoDerivative<double>& x);
FullPhysics::AutoDerivative<double> atan(const FullPhysics::AutoDerivative<double>& x);
}

%template(vector_auto_derivative) std::vector<FullPhysics::AutoDerivative<double> >;



