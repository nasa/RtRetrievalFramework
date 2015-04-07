// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "uq_ecmwf.h"
%}
%base_import(ecmwf);
%import "hdf_sounding_id.i"

%fp_shared_ptr(FullPhysics::UqEcmwf);

namespace FullPhysics {
class UqEcmwf : public Ecmwf  {
public:
    virtual ~UqEcmwf();
    UqEcmwf(const std::string& Fname);
};
}
