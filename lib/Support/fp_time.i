// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "fp_time.h"
%}
%base_import(generic_object)

//************************************************************
// Type map to use the Ruby type Time as input and output
//************************************************************

#ifdef SWIGRUBY
 %typemap(in) FullPhysics::Time {
  if(rb_obj_is_kind_of($input, rb_const_get(rb_cObject, 
					    rb_intern("Time"))) != Qtrue) {
    rb_raise(rb_eArgError, "Argument must be of type Time"); SWIG_fail;
  }
  $1 = FullPhysics::Time::time_unix(NUM2DBL(rb_funcall($input, rb_intern("to_f"), 0)));
 }

 %typemap(in) const FullPhysics::Time& (FullPhysics::Time t){
  if(rb_obj_is_kind_of($input, rb_const_get(rb_cObject, 
					    rb_intern("Time"))) != Qtrue) {
    rb_raise(rb_eArgError, "Argument must be of type Time"); SWIG_fail;
  }
  t = FullPhysics::Time::time_unix(NUM2DBL(rb_funcall($input, rb_intern("to_f"), 0)));
  $1 = &t;
 }

 %typemap(out) FullPhysics::Time {
  $result = rb_funcall(rb_const_get(rb_cObject, rb_intern("Time")), 
		       rb_intern("at"), 1, 
		       rb_float_new($1.unix_time()));
 }

 %typemap(out) const FullPhysics::Time& {
  $result = rb_funcall(rb_const_get(rb_cObject, rb_intern("Time")), 
		       rb_intern("at"), 1, 
		       rb_float_new($1->unix_time()));
 }

 %typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) FullPhysics::Time, const FullPhysics::Time& {
  $1 = (rb_obj_is_kind_of($input, rb_const_get(rb_cObject, 
					       rb_intern("Time"))) == Qtrue);
 }
#endif

// Don't think we actually want to use the Python typemap here, so
// we'll suppress this for now. Leave code here in case we change our
// mind.

//************************************************************
// Type map to use the Python type datetime as input and output. 
// We also support a python float as a unix time (seconds from
// January 1, 1970, as the normal unix time is specified).
//************************************************************

// #ifdef SWIGPYTHON

// %{
// //--------------------------------------------------------------
// // Function to return datetime module, importing if this is the
// // first call
// //--------------------------------------------------------------

// PyObject* datetime_module()
// {
//   static PyObject* mod = 0;
//   if(!mod)
//     mod = PyImport_ImportModule("datetime");
//   return mod;
// }

// //--------------------------------------------------------------
// // Return datetime.datetime object.
// //--------------------------------------------------------------

// PyObject* datetime_dot_datetime()
// {
//   static PyObject* obj = 0;
//   if(!obj)
//     obj = PyDict_GetItemString(PyModule_GetDict(datetime_module()) , 
// 			       "datetime");
//   return obj;
// }

// //--------------------------------------------------------------
// // Function to return time module, importing if this is the
// // first call
// //--------------------------------------------------------------

// PyObject* time_module()
// {
//   static PyObject* mod = 0;
//   if(!mod)
//     mod = PyImport_ImportModule("time");
//   return mod;
// }

// //--------------------------------------------------------------
// // Return time.mktime object.
// //--------------------------------------------------------------

// PyObject* time_dot_mktime()
// {
//   static PyObject* obj = 0;
//   if(!obj)
//     obj = PyDict_GetItemString(PyModule_GetDict(time_module()) , 
// 			       "mktime");
//   return obj;
// }

// //--------------------------------------------------------------
// // Convert datetime to FullPhysics::Time. This includes handling 
// // the microsecond field for fractional seconds (time.mktime
// // truncates to seconds).
// //--------------------------------------------------------------

// FullPhysics::Time datetime_to_time(PyObject* datetime)
// {
//   PyObject* tm = PyObject_CallMethod(datetime, 
// 				     const_cast<char*>("timetuple"), 
// 				     const_cast<char*>(""));
//   PythonObject res(PyObject_CallFunctionObjArgs(time_dot_mktime(), tm, NULL));
//   double tval = PyFloat_AsDouble(res);
//   PythonObject res2(PyObject_GetAttrString(datetime, "microsecond"));
//   tval += PyInt_AsLong(res2) / 1000000.0;
//   return FullPhysics::Time::time_unix(tval);
// }
// %}

// //--------------------------------------------------------------
// // We convert unix time to FullPhysics::Time object. 
// //--------------------------------------------------------------

// %typemap(in) FullPhysics::Time {
//   if(PyFloat_Check($input)) {
//     $1 = FullPhysics::Time::time_unix(PyFloat_AsDouble($input));
//   } else {
//     if(!PyObject_IsInstance($input, datetime_dot_datetime())) {
//       PyErr_SetString(PyExc_ValueError,"Expected a float or datetime object.");
//       return NULL;
//     }
//     $1 = datetime_to_time($input);
//   } 
// }

// //--------------------------------------------------------------
// // We convert unix time to FullPhysics::Time object.
// //--------------------------------------------------------------

// %typemap(in) const FullPhysics::Time& (FullPhysics::Time t){
//   if(PyFloat_Check($input)) {
//     t = FullPhysics::Time::time_unix(PyFloat_AsDouble($input));
//   } else {
//     if(!PyObject_IsInstance($input, datetime_dot_datetime())) {
//       PyErr_SetString(PyExc_ValueError,"Expected a float or datetime object.");
//       return NULL;
//     }
//     t = datetime_to_time($input);
//   } 
//   $1 = &t;
// }

// //--------------------------------------------------------------
// // Convert unix time from FullPhysics::Time object to python
// // datetime. 
// //--------------------------------------------------------------

// %typemap(out) FullPhysics::Time {
//   $result = PyObject_CallMethod(datetime_dot_datetime(), 
// 				const_cast<char*>("fromtimestamp"), 
// 				const_cast<char*>("d"), $1.unix_time());
// }

// //--------------------------------------------------------------
// // Convert unix time from FullPhysics::Time object to python
// // datetime. 
// //--------------------------------------------------------------

// %typemap(out) const FullPhysics::Time& {
//   $result = PyObject_CallMethod(datetime_dot_datetime(), 
// 				const_cast<char*>("fromtimestamp"), 
// 				const_cast<char*>("d"), $1->unix_time());
// }

// //--------------------------------------------------------------
// // Check if object is datetime or float.
// //--------------------------------------------------------------

// %typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) FullPhysics::Time, const FullPhysics::Time& {
//   $1 = PyObject_IsInstance($input, datetime_dot_datetime()) || 
//     PyFloat_Check($input);
// }
// #endif

%pythoncode {
import datetime
import time

def _new_time(pgs):
  return Time.time_pgs(pgs)
}

%fp_shared_ptr(FullPhysics::Time)

namespace FullPhysics {
class Time : public GenericObject {
public:
  static Time time_pgs(double pgs);
  static Time time_unix(double unix_time);
  %python_attribute(unix_time, double)
  %python_attribute(pgs_time, double)
  static Time parse_time(const std::string Time_string);
  std::string print_to_string() const;
  
#ifdef SWIGRUBY
  %extend {
// Ruby prefers defining this over defining operator<.
     int compare(const Time& T2) {
        if(*$self < T2)
         return -1;
        if(T2 < *$self)
         return 1;
        return 0;
     }
// Friendlier names for adding and subtracting
     Time plus(double T) { return *$self + T; }
     Time minus(double T) { return *$self - T; }
     double minus(const Time& T2) { return *$self - T2; }
  }
#endif
#ifdef SWIGPYTHON
  %extend {
     int __cmp__(const Time& T2) {
        if(*$self < T2)
         return -1;
        if(T2 < *$self)
         return 1;
        return 0;
     }
     Time __add__(double T) { return *$self + T; }
     Time __radd__(double T) { return *$self + T; }
     Time __sub__(double T) { return *$self - T; }
     double __sub__(const Time& T2) { return *$self - T2; }
  }
  %pythoncode {
def __reduce__(self):
  return _new_time, (self.pgs(),)

  }
#endif
};
}
