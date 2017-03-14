// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ecmwf.h"
%}
%base_import(meteorology)

%fp_shared_ptr(FullPhysics::Ecmwf);
%nodefaultctor FullPhysics::Ecmwf;

namespace FullPhysics {
class Ecmwf : public Meteorology {
public:
    virtual ~Ecmwf();
    %python_attribute(h2o_vmr, blitz::Array<double, 1>);
    %python_attribute(ozone_mmr, virtual blitz::Array<double, 1>)
    %python_attribute(ozone_vmr, blitz::Array<double, 1>)
    %python_attribute(windspeed_u, double);
    %python_attribute(windspeed_v, double);
    std::string print_to_string();
};
}

