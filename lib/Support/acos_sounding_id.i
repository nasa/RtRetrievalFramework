// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include <std_vector.i>
%include "common.i"

%{
#include "acos_sounding_id.h"
%}
%base_import(hdf_sounding_id);
%import "hdf_file.i"
%fp_shared_ptr(FullPhysics::AcosSoundingId);

namespace FullPhysics {
class AcosSoundingId : public HdfSoundingId {
public:
  enum SoundingType {S_SOUNDING, P_SOUNDING};
  virtual ~AcosSoundingId();
  AcosSoundingId(const HdfFile& File, const std::string& Sounding_id,
		SoundingType st);
  static std::vector<boost::shared_ptr<HdfSoundingId> >
  create(const HdfFile& File, const std::string& Sounding_id);
};
}

%template(vector_acos_sounding_id) std::vector<boost::shared_ptr<FullPhysics::HdfSoundingId> >;
