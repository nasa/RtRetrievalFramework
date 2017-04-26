%{
#include <boost/optional.hpp>
#include <blitz/array.h>
%}

%feature("novaluewrapper") boost::optional;

%typemap(out) boost::optional<blitz::Range> %{
    if (&$1) {
        $result = SWIG_NewPointerObj(new blitz::Range((&$1)->get()), $descriptor(blitz::Range*), SWIG_POINTER_OWN | 0);
    } else {
        $result = Py_None;
        Py_INCREF(Py_None);
    }
%}
