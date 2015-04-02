#ifndef UQ_SOUNDING_ID_H
#define UQ_SOUNDING_ID_H
#include "hdf_sounding_id.h"

namespace FullPhysics {
/****************************************************************//**
  This class reads a given file, and extracts out the sounding
  information. This determine the indexes into the HDF file for the
  given sounding (referred to as "frame number" and "sounding number").
*******************************************************************/
class UqSoundingId : public HdfSoundingId {
public:
  UqSoundingId(const HdfFile& File, const std::string& Sounding_id);
  virtual ~UqSoundingId() {}
  void print(std::ostream& Os) const;
private:
  void initialize(const HdfFile& File, const std::string& Sounding_id);
};

}
#endif
