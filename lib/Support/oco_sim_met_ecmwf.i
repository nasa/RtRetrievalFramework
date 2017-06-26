// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "oco_sim_met_ecmwf.h"
%}
%base_import(meteorology);
%import "hdf_sounding_id.i"

%fp_shared_ptr(FullPhysics::OcoSimMetEcmwf);

namespace FullPhysics {
class OcoSimMetEcmwf : public Meteorology {
public:
  virtual ~OcoSimMetEcmwf();
  OcoSimMetEcmwf(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& Hdf_sounding_id);
};
}
