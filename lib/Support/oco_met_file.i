// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "oco_met_file.h"
%}
%base_import(meteorology);
%import "hdf_sounding_id.i"

%fp_shared_ptr(FullPhysics::OcoMetFile);

namespace FullPhysics {
class OcoMetFile : public Meteorology  {
public:
  virtual ~OcoMetFile();
  OcoMetFile(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& Hdf_sounding_id);
};
}
