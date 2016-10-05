#ifndef OCO_SOUNDING_ID_H
#define OCO_SOUNDING_ID_H
#include "hdf_sounding_id.h"

namespace FullPhysics {
/****************************************************************//**
  This class reads a given file, and extracts out the sounding
  information. This determine the indexes into the HDF file for the
  given sounding (referred to as "frame number" and "sounding number").
*******************************************************************/
class OcoSoundingId : public HdfSoundingId {
public:
  OcoSoundingId(const HdfFile& File, const std::string& Sounding_id);
  /// The position in the OCO simulator scene file
  int scene_index() const
  { return number_sounding * frame_number() + sounding_number(); }
  virtual ~OcoSoundingId() {}
  void print(std::ostream& Os) const;
private:
  void initialize(const HdfFile& File, const std::string& Sounding_id);
  unsigned int sounding_pos_;
  int number_sounding;
};

}
#endif
