// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "hdf_sounding_id.h"
%}
%base_import(generic_object)

%fp_shared_ptr(FullPhysics::HdfSoundingId);
%nodefaultctor FullPhysics::HdfSoundingId;

namespace FullPhysics {
class HdfSoundingId : public GenericObject {
public:
  virtual ~HdfSoundingId();
  %python_attribute(frame_number, int);
  %python_attribute(sounding_number, int);
  %python_attribute(sounding_id, int64_t);
  std::string print_to_string() const;
};
}
