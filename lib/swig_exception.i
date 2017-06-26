// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

//--------------------------------------------------------------
// Translate exceptions into the appropriate language exception type
//--------------------------------------------------------------

%include "exception.i"

%exception {
  try {
    $action
  } catch (Swig::DirectorException &e) { 
    SWIG_fail; 
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%{
  // This is defined in swig_wrap.tmpl, so it gets put into swig_wrap.cc
  std::string parse_python_exception();
%}
%feature("director:except") {
    if ($error != NULL) {
      FullPhysics::Exception e;
      e << "Python error occured:\n"
	<< parse_python_exception();
      throw e;
    }
}

