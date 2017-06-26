// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

//--------------------------------------------------------------
// This provides common functions etc. used throughout our SWIG
// code. This should get included at the top of each swig file.
//--------------------------------------------------------------

// The module actually gets overridden by SWIG command line options
// when we build. But we still need to supply this to get the
// directors=1 and allprotected=1 set.

%module(directors="1", allprotected="1") full_physics

%{
#include <boost/shared_ptr.hpp>
#include <boost/rational.hpp>

//--------------------------------------------------------------
// Helper class for python that holds an object and when deleted
// decrements the reference to it.
//--------------------------------------------------------------

class PythonObject {
public:
  PythonObject(PyObject* Obj = 0) : obj(Obj) {}
  ~PythonObject() { Py_XDECREF(obj); }
  PyObject* obj;
  operator PyObject*() {return obj;}
};
%}

// Short cut for ingesting a base class
%define %base_import(NAME)
%import(module="full_physics_swig.NAME") "NAME.i"
%enddef

// Map std::string to and from the native string type
%naturalvar std::string;

// Include our own rules and common imports here.
%include "fp_shared_ptr.i"
%include "swig_exception.i"
%include "swig_print.i"
%include "swig_python_attribute.i"
%include "swig_pickle.i"
%import "swig_std.i"
%include "swig_array_inc.i"
%import "swig_array.i"
%import "swig_rational.i"

