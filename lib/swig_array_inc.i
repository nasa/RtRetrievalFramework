// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// Stuff to include everywhere for using blitz::array. We also have 
// swig_array.i which includes stuff compiled in one spot.

%{
// Don't want to use threads with ruby
#undef _REENTRANT
#include <blitz/array.h>
#include <blitz/range.h>
#define PY_ARRAY_UNIQUE_SYMBOL full_physics_ARRAY_API
#ifndef DO_IMPORT_ARRAY
#define NO_IMPORT_ARRAY
#endif
// See https://github.com/numpy/numpy/issues/3008 for explanation of
// this.
// We'll have to update this as the numpy API increases
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include "linear_algebra.h"
#include "fp_exception.h"

PyObject* numpy_module();
PyObject* numpy_dot_float64();
PyObject* numpy_dot_int32();
PyObject* numpy_dot_bool();

//--------------------------------------------------------------
// Helper routines to map a template type to the code numpy uses
// for that type.
//--------------------------------------------------------------

template<class T> int type_to_npy();
template<> inline int type_to_npy<double>() {return NPY_DOUBLE;}
template<> inline int type_to_npy<int>() {return NPY_INT;}
template<> inline int type_to_npy<bool>() {return NPY_BOOL;}

//--------------------------------------------------------------
// Use the numpy command "asarray" to convert various python 
// objects to a numpy object. This may return null, if the 
// "asarray" fails. 
//--------------------------------------------------------------

template<class T> PyObject* to_numpy(PyObject* obj);

template<> inline PyObject* to_numpy<double>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
					     PyString_FromString("asarray"), 
					     obj, numpy_dot_float64(), NULL);
  // Don't worry about errors , since we just return a null
  PyErr_Clear();
  return res;
}

template<> inline PyObject* to_numpy<bool>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
				    PyString_FromString("asarray"), 
				    obj, numpy_dot_bool(), NULL);
  PyErr_Clear();
  return res;
}

template<> inline PyObject* to_numpy<int>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
				    PyString_FromString("asarray"), 
				    obj, numpy_dot_int32(), NULL);
  PyErr_Clear();
  return res;
}

//--------------------------------------------------------------
// Convert a numpy array to a blitz::Array. The numpy should 
// already be the right data type before calling these (you can
// call to_numpy, if that is convenient). The underlying data is 
// still owned by the numpy object, so you need to make sure that
// the numpy object doesn't get deleted until you are done with
// the blitz::Array.
//
// If this fails, we throw an exception.
//--------------------------------------------------------------

template<class T, int D> inline blitz::Array<T, D> 
  to_blitz_array(PyObject* numpy_obj)
{
  PyArrayObject* numpy = (PyArrayObject*) numpy_obj;
  if(PyArray_NDIM(numpy) != D) {
    std::cerr << PyArray_NDIM(numpy) << "\n"
	      << D << "\n";
    throw 
      FullPhysics::Exception("Dimension of array is not the expected size");
  }
  if(PyArray_TYPE(numpy) != type_to_npy<T>()) {
    throw 
      FullPhysics::Exception("Type of array not the expected type");
  }
  blitz::TinyVector<int, D> shape, stride;
  for(int i = 0; i < D; ++i) {
    shape(i) = PyArray_DIM(numpy, i);
    // Note numpy stride is in terms of bytes, while blitz in in terms
    // of type T.
    stride(i) = PyArray_STRIDE(numpy, i) / sizeof(T);
    if((int) (stride(i) * sizeof(T)) != (int) PyArray_STRIDE(numpy, i)) {
      throw 
	FullPhysics::Exception("blitz::Array can't handle strides that aren't an even multiple of sizeof(T)");
    }
  }
  return blitz::Array<T, D>((T*)PyArray_DATA(numpy), shape, stride, 
			    blitz::neverDeleteData);
}

%}

