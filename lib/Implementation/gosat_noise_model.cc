#include "gosat_noise_model.h"

#include <gsl/gsl_interp.h>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(GosatNoiseModel, NoiseModel)
.def(luabind::constructor<const HeritageFile&>())
.def(luabind::constructor<const HeritageFile&,const HeritageFile&>())
.def(luabind::constructor<const HdfFile&,const HdfSoundingId&, const Instrument&>())
.def(luabind::constructor<const HdfFile&,const HdfSoundingId&, const std::vector<std::string>&>())
.def(luabind::constructor<const HdfFile&,const HdfSoundingId&, const Instrument&,const HeritageFile&>())
.def(luabind::constructor<const HdfFile&,const HdfSoundingId&, const std::vector<std::string>&,const HeritageFile&>())
.def(luabind::constructor<const HdfFile&,const HdfSoundingId&, const Instrument&,const HdfFile&>())
.def(luabind::constructor<const HdfFile&,const HdfSoundingId&, const std::vector<std::string>&,const HdfFile&>())
.def(luabind::constructor<const HdfFile&,const HdfSoundingId&, const Instrument&,const HdfFile&,const std::string&>())
.def(luabind::constructor<const HdfFile&,const HdfSoundingId&, const std::vector<std::string>&,const HdfFile&,const std::string&>())
REGISTER_LUA_END()
#endif


//-----------------------------------------------------------------------
/// Creates a new GOSAT noise model from ascii heritage file.
//-----------------------------------------------------------------------

GosatNoiseModel::GosatNoiseModel(const HeritageFile& Noise_file) {
  read_ascii_noise_file(Noise_file);
}

//-----------------------------------------------------------------------
/// Creates a new GOSAT noise model from ascii heritage file and uses
/// empirical noise coefficients from the supplied file
//-----------------------------------------------------------------------

GosatNoiseModel::GosatNoiseModel(const HeritageFile& Noise_file, const HeritageFile& Emp_Coeff_File) {
  read_ascii_noise_file(Noise_file);
  read_empirical_noise_coeffs(Emp_Coeff_File);
}

//-----------------------------------------------------------------------
/// Creates a new GOSAT noise model based on contents of HDF file data
/// at the given exposure index
//-----------------------------------------------------------------------

GosatNoiseModel::GosatNoiseModel(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, const Instrument& Inst) {
  std::vector<std::string> Hdf_band_name;
  for(int spec_idx = 0; spec_idx < Inst.number_spectrometer(); spec_idx++)
    Hdf_band_name.push_back(Inst.hdf_band_name(spec_idx));
  read_hdf_noise_file(Noise_file, Sounding_Id, Hdf_band_name);
}

//-----------------------------------------------------------------------
/// Creates a new GOSAT noise model based on contents of HDF file data
/// at the given exposure index
//-----------------------------------------------------------------------
GosatNoiseModel::GosatNoiseModel(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, const std::vector<std::string>& Hdf_band_name) {
  read_hdf_noise_file(Noise_file, Sounding_Id, Hdf_band_name);
}

//-----------------------------------------------------------------------
/// Creates a new GOSAT noise model based on contents of HDF file data
/// with empirical noise correction values read in
//-----------------------------------------------------------------------

GosatNoiseModel::GosatNoiseModel(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, const Instrument& Inst, const HeritageFile& Emp_Coeff_File) {
  std::vector<std::string> Hdf_band_name;
  for(int spec_idx = 0; spec_idx < Inst.number_spectrometer(); spec_idx++)
    Hdf_band_name.push_back(Inst.hdf_band_name(spec_idx));
  read_hdf_noise_file(Noise_file, Sounding_Id, Hdf_band_name);
  read_empirical_noise_coeffs(Emp_Coeff_File);
}

//-----------------------------------------------------------------------
/// Creates a new GOSAT noise model based on contents of HDF file data
/// with empirical noise correction values read in
//-----------------------------------------------------------------------

GosatNoiseModel::GosatNoiseModel(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, const std::vector<std::string>& Hdf_band_name, const HeritageFile& Emp_Coeff_File) {
  read_hdf_noise_file(Noise_file, Sounding_Id, Hdf_band_name);
  read_empirical_noise_coeffs(Emp_Coeff_File);
}

//-----------------------------------------------------------------------
/// Creates a new GOSAT noise model based on contents of HDF file data
/// with empirical noise correction values read in from the HDF config
/// file.
//-----------------------------------------------------------------------

GosatNoiseModel::GosatNoiseModel
(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, 
 const Instrument& Inst, const HdfFile& Emp_Coeff_File, 
 const std::string& Group_name) {
  std::vector<std::string> Hdf_band_name;
  for(int spec_idx = 0; spec_idx < Inst.number_spectrometer(); spec_idx++)
    Hdf_band_name.push_back(Inst.hdf_band_name(spec_idx));
  read_hdf_noise_file(Noise_file, Sounding_Id, Hdf_band_name);
  read_empirical_noise_coeffs(Emp_Coeff_File, Group_name);
}

//-----------------------------------------------------------------------
/// Creates a new GOSAT noise model based on contents of HDF file data
/// with empirical noise correction values read in from the HDF config
/// file.
//-----------------------------------------------------------------------

GosatNoiseModel::GosatNoiseModel
(const HdfFile& Noise_file, const HdfSoundingId& Sounding_Id, 
 const std::vector<std::string>& Hdf_band_name, const HdfFile& Emp_Coeff_File, 
 const std::string& Group_name) {
  read_hdf_noise_file(Noise_file, Sounding_Id, Hdf_band_name);
  read_empirical_noise_coeffs(Emp_Coeff_File, Group_name);
}

//-----------------------------------------------------------------------

void GosatNoiseModel::read_hdf_noise_file(const HdfFile& Noise_file, 
  const HdfSoundingId& Sounding_Id, const std::vector<std::string>& Hdf_band_name) {
 
  firstIndex i1; 
  secondIndex i2;
  thirdIndex i3; 

  int frame_idx = Sounding_Id.frame_number();
  int snd_idx = Sounding_Id.sounding_number();

  Array<std::string, 2> ds_gain_codes = 
    Noise_file.read_field<std::string, 2>("SoundingHeader/gain_swir", 
					  TinyVector<int, 2>(frame_idx, snd_idx), 
					  TinyVector<int, 2>(1,1));
  char gain_code = ds_gain_codes(0, 0)[0];

  // Determine the prefix of the dataset name used for extracting CNV data
  // based on L1B file's gain code
  string cnv_ds_name_beg;
  switch(gain_code) {
  case 'H':
    cnv_ds_name_beg = "InstrumentHeader/cnv_coef_highgain_";
    break;
  case 'M':
    cnv_ds_name_beg = "InstrumentHeader/cnv_coef_medgain_";
    break;
  default:
    Exception e;
    e << "Gain code: " << gain_code << " for sounding id: " << Sounding_Id.sounding_id() << " not supported.";
    throw e;
  }

  // Determine the maximum size needed for CNV data
  int max_cnv_color_size = 0;
  for(int spec_idx = 0; spec_idx < (int) Hdf_band_name.size(); spec_idx++) {
    stringstream ds_name;
    ds_name << cnv_ds_name_beg << Hdf_band_name[spec_idx];
    max_cnv_color_size = max(max_cnv_color_size, Noise_file.read_shape<3>(ds_name.str())[thirdDim]);
  }

  // Now read per aband data, also set up inp_color_index to just
  // contain the indexes of valid locations for the data
  band_noise_val.resize(Hdf_band_name.size());
  inp_color_index.resize(max_cnv_color_size, Hdf_band_name.size());
  coef_rad_cnv.resize(max_cnv_color_size, Hdf_band_name.size());
 
  inp_color_index = 0; // Reset all to zero since this one will be searched for indexes
  for(int spec_idx = 0; spec_idx < (int) Hdf_band_name.size(); spec_idx++) {
    stringstream ds_name;
    ds_name << cnv_ds_name_beg << Hdf_band_name[spec_idx];
    
    int file_cnv_data_size = Noise_file.read_shape<3>(ds_name.str())[thirdDim];
    Range cnv_range = Range(0,file_cnv_data_size-1);
    coef_rad_cnv(cnv_range, spec_idx) = Noise_file.read_field<double, 3>
      (ds_name.str(), 
       TinyVector<int, 3>(frame_idx, snd_idx, 0), 
       TinyVector<int, 3>(1,1,file_cnv_data_size))(0,0,Range::all());

    for (int color_idx = 0; color_idx < file_cnv_data_size; 
	 color_idx++) {
      inp_color_index(color_idx, spec_idx) = color_idx;
    }

    // Build 2.9.00 of L1B products from ACOS renamed noise_XXX dataset to noise_XXX_l1b
    // Check for this falling back on the pre 2.9.00 dataset name
    try {
      band_noise_val(spec_idx) =
	Noise_file.read_field<double, 2>("SoundingSpectra/noise_" + Hdf_band_name[spec_idx] + "_l1b", 
					 TinyVector<int, 2>(frame_idx, snd_idx),
					 TinyVector<int, 2>(1, 1))(0, 0);
    } catch(Exception not_found_error ) {
      band_noise_val(spec_idx) =
	Noise_file.read_field<double, 2>("SoundingSpectra/noise_" + Hdf_band_name[spec_idx], 
					 TinyVector<int, 2>(frame_idx, snd_idx),
					 TinyVector<int, 2>(1, 1))(0, 0);
    }
  }
}

void GosatNoiseModel::read_ascii_noise_file(const HeritageFile& Noise_file) {
  Range all = Range::all();
  
  std::vector<double> noise_vals_vec = Noise_file.value<std::vector<double> >("band_noise_values");
  band_noise_val.resize(noise_vals_vec.size());
  for(unsigned int i = 0; i < noise_vals_vec.size(); i++)
    band_noise_val(i) = noise_vals_vec[i];

  int max_coef_vals = Noise_file.data().rows();
  
  inp_color_index.resize(max_coef_vals, band_noise_val.rows());
  coef_rad_cnv.resize(max_coef_vals, band_noise_val.rows());

  for(int band_idx = 0; band_idx < band_noise_val.rows(); band_idx++) {
    std::stringstream color_col_name, cnv_col_name;
    color_col_name << "COLOR_" << band_idx+1;
    cnv_col_name << "CNV_COEF_" << band_idx+1;

    // Convert from fortran to C++ index with the -1 here
    inp_color_index(all, band_idx) = Noise_file.data(color_col_name.str()) - 1;
    coef_rad_cnv(all, band_idx) = Noise_file.data(cnv_col_name.str());
  }
}

void GosatNoiseModel::read_empirical_noise_coeffs(const HeritageFile& Emp_Coeff_File) {

  // Dimensions should be: coeff x spectrometer
  empirical_noise_coef.reference(Emp_Coeff_File.data().copy());
}

void GosatNoiseModel::read_empirical_noise_coeffs
(const HdfFile& Emp_Coeff_File, const std::string& Group_name) 
{
  empirical_noise_coef.reference
    (Emp_Coeff_File.read_field<double, 2>(Group_name + "/Empirical_Noise"));
}



//-----------------------------------------------------------------------
/// Calculate uncertainty related to GOSAT radiance using the
/// RadCnv values and noise value associated with a certain
/// sounding.
//-----------------------------------------------------------------------

blitz::Array<double, 1> GosatNoiseModel::uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const {
  
  Array<double, 1> uncert_res(Radiance.shape());

  // Only send the indexes for this spectrometer, of the length of the radiance to
  // the gsl search routine
  Array<double, 1> spec_color(Radiance.rows());
  spec_color = inp_color_index(Range(0, Radiance.rows()-1), Spec_index);

  // Compute signal level
  double signal_level = 0;
  if(empirical_noise_coef.rows() > 0) {
    /// \todo Move signal computation to a shared location with 
    /// the computation also used by error analysis
    const int nrad = 10;
    Array<double, 1> rad(Radiance.rows());
    rad = Radiance(Range::all());
    std::sort(rad.data(), rad.data() + Radiance.rows()); // Min to max value
    rad.reverseSelf(firstDim);	     // Now max to min value
    Range r2(0, std::min(nrad - 1, rad.rows() - 1));
    signal_level = sum(rad(r2) / r2.length());
  }
  
  for (int pix_idx = 0; pix_idx < Radiance.rows(); pix_idx++) {
    int color_index = gsl_interp_bsearch(spec_color.dataFirst(), pix_idx, 0, spec_color.rows());
    uncert_res(pix_idx) = coef_rad_cnv(color_index, Spec_index) * band_noise_val(Spec_index);

    // Apply empircal noise coefficients if loaded
    if(empirical_noise_coef.rows() > 0) {
      uncert_res(pix_idx) = empirical_noise_coef(0, Spec_index) * signal_level + \
	empirical_noise_coef(1, Spec_index) * uncert_res(pix_idx);
    }
  }

  return uncert_res;
}

void GosatNoiseModel::print(std::ostream& Os) const
{
  Os << "Gosat Noise Model:" << std::endl
     << "  Num Bands:    " << coef_rad_cnv.cols() << std::endl
     << "  Max Colors:   " << coef_rad_cnv.rows() << std::endl
     << "  Noise Values: " << band_noise_val << std::endl;
}
