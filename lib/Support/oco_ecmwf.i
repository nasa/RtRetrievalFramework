// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "oco_ecmwf.h"
%}
%base_import(meteorology);
%import "hdf_sounding_id.i"

%fp_shared_ptr(FullPhysics::OcoEcmwf);

namespace FullPhysics {
class OcoEcmwf : public Meteorology  {
public:
  virtual ~OcoEcmwf();
  OcoEcmwf(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& Hdf_sounding_id);
};
}
