#ifndef HDF_SOUNDING_ID_H
#define HDF_SOUNDING_ID_H
#include "hdf_file.h"
#include "heritage_file.h"
#include "printable.h"
#include <stdint.h>

namespace FullPhysics {
/****************************************************************//**
  This class defines the interface for a sounding id referring
  to data within an HDF file. This amounts to translating
  a string into frame and sounding indexes where data should
  be read from the HDF files. The details of the exact meaning
  of the string and how the indexing is computed may differ
  from file type to file type.
*******************************************************************/
class HdfSoundingId : public Printable<HdfSoundingId> {
public:
    virtual ~HdfSoundingId() {}

    virtual int frame_number() const {return frame_number_;}
    virtual int sounding_number() const {return sounding_number_;}
    virtual int64_t sounding_id() const { return sounding_id_int; }
    virtual void print(std::ostream& Os) const = 0;
protected:
    int frame_number_;
    int sounding_number_;
    int64_t sounding_id_int;
};

}
#endif
