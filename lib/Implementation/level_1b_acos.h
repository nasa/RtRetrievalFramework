#ifndef LEVEL_1B_ACOS_H
#define LEVEL_1B_ACOS_H
#include "level_1b_hdf.h"
#include "acos_sounding_id.h"
#include "fp_exception.h"
#include "hdf_file.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This reads a Level 1B file that is in the HDF format.
*******************************************************************/
class Level1bAcos: public Level1bHdf {
public:
  Level1bAcos(const std::string& Fname, 
	      const boost::shared_ptr<HdfSoundingId>& Sounding_id);
  Level1bAcos(const boost::shared_ptr<HdfFile>& Hfile, 
	      const boost::shared_ptr<HdfSoundingId>& Sounding_id);

//-----------------------------------------------------------------------
/// Percentage of land in sounding Field of View
//-----------------------------------------------------------------------

  double land_fraction(int spec_index) const {
    range_check(spec_index, 0, number_spectrometer());
    return land_fraction_(spec_index); 
  }

//-----------------------------------------------------------------------
/// True if this is H gain
//-----------------------------------------------------------------------

  bool is_h_gain() const { return is_h_gain_; }

//-----------------------------------------------------------------------
/// True if this is M gain
//-----------------------------------------------------------------------

  bool is_m_gain() const { return !is_h_gain(); }
protected:
  virtual SpectralRange radiance_no_uncertainty(int Spec_index) const;
private:
  bool is_h_gain_;
  blitz::Array<double,1> land_fraction_;
  void initialize();
};
}
#endif
