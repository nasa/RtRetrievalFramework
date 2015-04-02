#ifndef OCO_NOISE_MODEL_H
#define OCO_NOISE_MODEL_H
#include "hdf_file.h"
#include "hdf_sounding_id.h"
#include "noise_model.h"

namespace FullPhysics {
/****************************************************************//**
  This class creates a OCO noise model using inputs from the supplied file
*******************************************************************/
class OcoNoiseModel: public NoiseModel {
public:

  /// Creates a new OCO noise model from the L1B HDF5 file
  OcoNoiseModel(const HdfFile& Hfile, 
		const HdfSoundingId& Sounding_id) 
  { max_ms_.resize(3);
    max_ms_ = 7.0e20, 2.45e20, 1.25e20;
    read_hdf_noise(Hfile, Sounding_id); }

  OcoNoiseModel(const HdfFile& Hfile, 
		const HdfSoundingId& Sounding_id,
		const blitz::Array<double, 1>& Max_meas_signal) 
    : max_ms_(Max_meas_signal)
  { read_hdf_noise(Hfile, Sounding_id); }

  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const;
    
  virtual void print(std::ostream& Os) const;
  
  blitz::Array<double, 1> coef_photon(int Spec_index) 
  { range_max_check(Spec_index, coef_photon_.extent(blitz::secondDim)+1);
    return coef_photon_(blitz::Range::all(), Spec_index); }
  
  blitz::Array<double, 1> coef_background(int Spec_index) 
  { range_max_check(Spec_index, coef_background_.extent(blitz::secondDim)+1);
    return coef_background_(blitz::Range::all(), Spec_index); }

protected:
  // For inheriting class
  OcoNoiseModel() {}

  blitz::Array<double, 2> coef_photon_;
  blitz::Array<double, 2> coef_background_;
  blitz::Array<double, 1> max_ms_;

private:
  virtual void read_hdf_noise(const HdfFile& Hfile, const HdfSoundingId& Sounding_id);
 
};
}
#endif


