// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "acos_met_file.h"
%}
%base_import(meteorology);
%import "hdf_sounding_id.i"

%fp_shared_ptr(FullPhysics::AcosMetFile);

namespace FullPhysics {
class AcosMetFile : public Meteorology  {
public:
  virtual ~AcosMetFile();
  AcosMetFile(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& 
	    Hdf_sounding_id, bool Avg_sounding_number);
};
}
