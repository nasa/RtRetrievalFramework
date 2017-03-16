// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "uq_ecmwf.h"
%}
%base_import(meteorology);
%import "hdf_sounding_id.i"

%fp_shared_ptr(FullPhysics::UqEcmwf);

namespace FullPhysics {
class UqEcmwf : public Meteorology {
public:
    virtual ~UqEcmwf();
    UqEcmwf(const std::string& Fname);
};
}
