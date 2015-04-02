// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "uq_sounding_id.h"
%}
%base_import(hdf_sounding_id);
%import "hdf_file.i"

%fp_shared_ptr(FullPhysics::UqSoundingId);

namespace FullPhysics {
class UqSoundingId : public HdfSoundingId {
public:
    UqSoundingId(const HdfFile& File, const std::string& Sounding_id);
    virtual ~UqSoundingId();
};
}
