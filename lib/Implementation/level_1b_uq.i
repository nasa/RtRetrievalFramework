// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "level_1b_uq.h"
%}
%base_import(level_1b_oco)
%import "hdf_sounding_id.i"
%import "hdf_file.i"
%fp_shared_ptr(FullPhysics::Level1bUq);

namespace FullPhysics {
class Level1bUq: public Level1bOco {
public:
    Level1bUq(const boost::shared_ptr<HdfFile>& Hfile, 
              const boost::shared_ptr<HdfSoundingId>& Sounding_id);
    void set_radiance(int Spec_index, boost::shared_ptr<SpectralRange>& Rad);
};
}
