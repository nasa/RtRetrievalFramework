// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "acos_ecmwf.h"
%}
%base_import(meteorology);
%import "hdf_sounding_id.i"

%fp_shared_ptr(FullPhysics::AcosEcmwf);

namespace FullPhysics {
class AcosEcmwf : public Meteorology  {
public:
  virtual ~AcosEcmwf();
  AcosEcmwf(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& 
	    Hdf_sounding_id, bool Avg_sounding_number);
};
}
