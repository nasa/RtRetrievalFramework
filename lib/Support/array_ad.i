// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "array_ad.h"
%}

%import "auto_derivative.i"

namespace FullPhysics {
template<class T, int D> class ArrayAd;
}

%pythoncode %{
import numpy as np

def np_to_array_ad(a):
    '''Convert a numpy array of AutoDerivatives to a ArrayAd'''
    nvar = 0
    for i in a.flat:
        if(i.number_variable > 0):
            nvar = i.number_variable
            break
    if(len(a.shape) == 1):
        res = ArrayAd_double_1(a.shape[0], nvar)
        for i1 in range(res.rows):
            res[i1]  = a[i1]
    elif(len(a.shape) == 2):
        res = ArrayAd_double_2(a.shape[0], a.shape[1], nvar)
        for i1 in range(res.rows):
            for i2 in range(res.cols):
                res[i1,i2]  = a[i1, i2]
    elif(len(a.shape) == 3):
        res = ArrayAd_double_3(a.shape[0], a.shape[1], a.shape[2], nvar)
        for i1 in range(res.rows):
            for i2 in range(res.cols):
                for i3 in range(res.depth):
                    res[i1,i2, i3]  = a[i1, i2, i3]
    else:
        raise IndexError("np_to_array_ad only implemented for dimension <= 3")
    return res
%}
%define %array_ad_template(NAME, TYPE, DIM, DIMP1)

namespace FullPhysics {
template<> class ArrayAd<TYPE, DIM>
{
public:
  ArrayAd<TYPE, DIM>(const blitz::Array<AutoDerivative<TYPE>, DIM>& V);
  ArrayAd<TYPE, DIM>(int n1, int nvar);
  ArrayAd<TYPE, DIM>(int n1, int n2, int nvar);
  ArrayAd<TYPE, DIM>(int n1, int n2, int n3, int nvar);
  ArrayAd<TYPE, DIM>(int n1, int n2, int n3, int n4, int nvar);
  ArrayAd<TYPE, DIM>(int n1, int n2, int n3, int n4, int n5, int nvar);
  ArrayAd<TYPE, DIM>(const blitz::Array<TYPE, DIM>& FORCE_COPY, 
		     const blitz::Array<TYPE, DIMP1>& FORCE_COPY,
		     bool Is_const = false);
  ArrayAd<TYPE, DIM>(const blitz::Array<TYPE, DIM>& FORCE_COPY);
  void resize_number_variable(int nvar);
  void resize(int n1, int nvar);
  void resize(int n1, int n2, int nvar);
  void resize(int n1, int n2, int n3, int nvar);
  void resize(int n1, int n2, int n3, int n4, int nvar);
  void resize(int n1, int n2, int n3, int n4, int n5, int nvar);
  %python_attribute(value, blitz::Array<TYPE, DIM>)
  %python_attribute(jacobian, blitz::Array<TYPE, DIMP1>)
  %python_attribute(rows, int)
  %python_attribute(cols, int)
  %python_attribute(depth, int)
  %python_attribute(is_constant, bool)
  %python_attribute(number_variable, int)
  void reference(const FullPhysics::ArrayAd<TYPE, DIM>& V);
  FullPhysics::ArrayAd<TYPE, DIM> copy() const;
  // When we pass multiple items in, this gets set to one tuple
  // argument. This confuses SWIG, so add one level of indirection
  // where we unpack this.
%pythoncode %{
def __array__(self):
    if(DIM == 1):
        res = np.empty([self.rows], np.object)
        for i1 in range(self.rows):
            res[i1] = self[i1]
    elif(DIM ==2):
        res = np.empty([self.rows, self.cols], np.object)
        for i1 in range(self.rows):
            for i2 in range(self.cols):
                res[i1,i2] = self[i1, i2]
    elif(DIM ==3):
        res = np.empty([self.rows, self.cols, self.depth], np.object)
        for i1 in range(self.rows):
            for i2 in range(self.cols):
                for i3 in range(self.depth):
                    res[i1,i2, i3] = self[i1, i2, i3]
    else:
      raise IndexError("__array__ only implemented to dimensions <= 3")
    return res

def slice_data(self, key):
  if not type(key) is tuple:
    key = (key,)

  ad_val = self.value[key]
  ad_jac = self.jacobian[key + (slice(None),)]

  num_dim = len(ad_val.shape)

  return eval("ArrayAd_double_%d" % num_dim)(ad_val, ad_jac)

def __getitem__(self, index):
  if(DIM == 1):
    if type(index) is slice:
      return self.slice_data(index)
    else:
      return self.read(index)
  else:
    if any(type(x) is slice for x in index):
      return self.slice_data(index)
    else:
      return self.read(*index)

def __setitem__(self, index, val):
  if(DIM == 1):
    if type(index) is slice:
      raise NotImplementedError("__setitem__ can not be used for setting values to slices")
    self.write(index, val)
  else:
    if any(type(x) is slice for x in index):
      raise NotImplementedError("__setitem__ can not be used for setting values to slices")
    t = list(index)
    t.append(val)
    self.write(*t)
%}
  %extend {
    AutoDerivative<TYPE> read(int i1) const { return (*$self)(i1); }
    AutoDerivative<TYPE> read(int i1, int i2) const { return (*$self)(i1, i2); }
    AutoDerivative<TYPE> read(int i1, int i2, int i3) const { return (*$self)(i1, i2, i3); }
    AutoDerivative<TYPE> read(int i1, int i2, int i3, int i4) const 
    { return (*$self)(i1, i2, i3, i4); }

    void write(int i1, const AutoDerivative<TYPE>& V) 
    { (*$self)(i1) = V; }
    void write(int i1, int i2, const AutoDerivative<TYPE>& V) 
    { (*$self)(i1, i2) = V; }
    void write(int i1, int i2, int i3, const AutoDerivative<TYPE>& V) 
    { (*$self)(i1, i2, i3) = V; }
    void write(int i1, int i2, int i3, int i4, const AutoDerivative<TYPE>& V)
    { (*$self)(i1, i2, i3, i4) = V; }
  }
  // 1 here is the pickle format version, so we can tell if we try to
  // read data with a different format version than the code here.
  %pickle_init(1, self.value, self.jacobian, self.is_constant)
 };
 
}

%template(NAME) FullPhysics::ArrayAd<TYPE, DIM>;
%enddef
%array_ad_template(ArrayAd_double_1, double, 1, 2);
%array_ad_template(ArrayAd_double_2, double, 2, 3);
%array_ad_template(ArrayAd_double_3, double, 3, 4);
%array_ad_template(ArrayAd_double_4, double, 4, 5);


