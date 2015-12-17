// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "rayleigh_greek_moment.h"
%}

namespace FullPhysics {
class RayleighGreekMoment {
public:
  static blitz::Array<double, 2> array();
};
}
