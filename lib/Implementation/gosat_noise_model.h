#ifndef GOSAT_NOISE_MODEL_H
#define GOSAT_NOISE_MODEL_H
#include "heritage_file.h"
#include "hdf_file.h"
#include "hdf_sounding_id.h"
#include "instrument.h"

#include "noise_model.h"

namespace FullPhysics {
/****************************************************************//**
  This class creates a GOSAT noise model using inputs from the supplied file
  which can be either a hertitage file format or HDF file
*******************************************************************/
class GosatNoiseModel: public NoiseModel {
public:

  GosatNoiseModel(const HeritageFile& Noise_file);
  GosatNoiseModel(const HeritageFile& Noise_file, const HeritageFile& Emp_Coeff_File);

  GosatNoiseModel(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, 
		  const Instrument& Inst);
  GosatNoiseModel(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, 
		  const std::vector<std::string>& Hdf_band_name);

  GosatNoiseModel(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, 
		  const Instrument& Inst, const HeritageFile& Emp_Coeff_File);
  GosatNoiseModel(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, 
		  const std::vector<std::string>& Hdf_band_name, const HeritageFile& Emp_Coeff_File);

  GosatNoiseModel(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, 
		  const Instrument& Inst, const HdfFile& Emp_Coeff_File,
		  const std::string& Group_name = "Instrument");
  GosatNoiseModel(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, 
		  const std::vector<std::string>& Hdf_band_name, const HdfFile& Emp_Coeff_File,
		  const std::string& Group_name = "Instrument");

  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const;

  virtual void print(std::ostream& Os) const;

private:

  void read_hdf_noise_file(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, const std::vector<std::string>& Hdf_band_name);
  void read_ascii_noise_file(const HeritageFile& Noise_file);
  void read_empirical_noise_coeffs(const HeritageFile& Emp_Coeff_File);
  void read_empirical_noise_coeffs(const HdfFile& Emp_Coeff_File, 
				   const std::string& Group_name);

  blitz::Array<double, 2> empirical_noise_coef;

  blitz::Array<double, 1> band_noise_val;
  blitz::Array<double, 2> inp_color_index;
  blitz::Array<double, 2> coef_rad_cnv;
  
};
}
#endif


