// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "closest_point.h"
%}

namespace FullPhysics {

template <class T1, class T2> 
int closest_point(const blitz::Array<T1, 1>& target, T2 data);

template <class T1, class T2> 
blitz::Array<int, 1> closest_point(const blitz::Array<T1, 1>& target,
                                   const blitz::Array<T2, 1>& data);

}
