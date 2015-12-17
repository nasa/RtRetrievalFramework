#ifndef ACOS_SOUNDING_ID_H
#define ACOS_SOUNDING_ID_H
#include "hdf_sounding_id.h"

namespace FullPhysics {
/****************************************************************//**
  This class reads a given file, and extracts out the sounding
  information. This determine the index into the HDF file for the
  given sounding (referred to as "frame number"), and determines which
  sounding we are using (s or p), or if we are doing averaging.
*******************************************************************/
class AcosSoundingId : public HdfSoundingId {
public:
  enum SoundingType {S_SOUNDING, P_SOUNDING};
  virtual ~AcosSoundingId() {}
  AcosSoundingId(const HdfFile& File, const std::string& Sounding_id,
		 SoundingType Sounding_type);
  static std::vector<boost::shared_ptr<HdfSoundingId> >
  create(const HdfFile& File, const std::string& Sounding_id);
  void print(std::ostream& Os) const;
private:
  void initialize(const HdfFile& File, const std::string& Sounding_id,
		  SoundingType Sounding_type);
};

}
#endif
