// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "level_1b_acos.h"
%}
%base_import(level_1b_hdf)
%import "hdf_sounding_id.i"
%import "hdf_file.i"
%import "spectral_range.i"
%fp_shared_ptr(FullPhysics::Level1bAcos);

namespace FullPhysics {
class Level1bAcos: public Level1bHdf {
public:
  Level1bAcos(const std::string& Fname, 
	      const boost::shared_ptr<HdfSoundingId>& Sounding_id);
  Level1bAcos(const boost::shared_ptr<HdfFile>& Hfile, 
	      const boost::shared_ptr<HdfSoundingId>& Sounding_id);
  virtual SpectralRange radiance(int Spec_index) const;
  double land_fraction(int spec_index) const;
  bool is_h_gain() const;
  bool is_m_gain() const;
};
}
