//--------------------------------------------------------------
// Function that allows us to access something as a python 
// attribute.
//--------------------------------------------------------------

 // We stick the TYPE at the end with the "..." so that we can pass
 // templates that contain "," through. This avoid the normal problem
 // with the preprocessor of getting confused by "," in a template
 // name, thinking these are separate arguments to macro

%define %python_attribute(NAME, TYPE...)
  %rename(_v_ ## NAME) NAME;
  TYPE NAME() const;
%pythoncode {
@property
def NAME(self):
    return self._v_ ## NAME()
}
%enddef

%define %python_attribute_with_set(NAME, TYPE...)
  %rename(_v_ ## NAME) NAME;
  TYPE NAME() const;
  void NAME(const TYPE& V);
%pythoncode {
@property
def NAME(self):
    return self._v_ ## NAME()

@NAME.setter
def NAME(self, value):
  self._v_ ## NAME(value)
}
%enddef

%define %python_attribute_with_set_ptr(NAME, TYPE...)
  %rename(_v_ ## NAME) NAME;
  boost::shared_ptr<TYPE> NAME() const;
  void NAME(const TYPE& V);
%pythoncode {
@property
def NAME(self):
    return self._v_ ## NAME()

@NAME.setter
def NAME(self, value):
  self._v_ ## NAME(value)
}
%enddef

%define %python_attribute_nonconst(NAME, TYPE...)
  %rename(_v_ ## NAME) NAME;
  TYPE NAME();
%pythoncode {
@property
def NAME(self):
    return self._v_ ## NAME()
}
%enddef

%define %python_attribute2(NAME, NAME2, TYPE...)
  %rename(_v_ ## NAME) NAME2;
  TYPE NAME2() const;
%pythoncode {
@property
def NAME(self):
    return self._v_ ## NAME()
}
%enddef

%define %python_attribute2_nonconst(NAME, NAME2, TYPE...)
  %rename(_v_ ## NAME) NAME2;
  TYPE NAME2();
%pythoncode {
@property
def NAME(self):
    return self._v_ ## NAME()
}
%enddef

%define %python_attribute_abstract(NAME, TYPE...)
  %rename(_v_ ## NAME) NAME;
  virtual TYPE NAME() const = 0;
%pythoncode {
@property
def NAME(self):
    return self._v_ ## NAME()
}

%enddef

%define %python_attribute_derived(NAME, TYPE...)
  %rename(_v_ ## NAME) NAME;
  virtual TYPE NAME() const;
%enddef

