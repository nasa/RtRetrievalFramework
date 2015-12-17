// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "oco_sounding_id.h"
%}
%base_import(hdf_sounding_id);
%import "hdf_file.i"

%fp_shared_ptr(FullPhysics::OcoSoundingId);

namespace FullPhysics {
class OcoSoundingId : public HdfSoundingId {
public:
  OcoSoundingId(const HdfFile& File, const std::string& Sounding_id);
  virtual ~OcoSoundingId();
};
}
