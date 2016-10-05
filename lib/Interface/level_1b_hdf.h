#ifndef LEVEL_1B_HDF_H
#define LEVEL_1B_HDF_H
#include "level_1b.h"
#include "hdf_sounding_id.h"
#include "fp_exception.h"
#include "hdf_file.h"
#include "noise_model.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This is an abstract interface class which mainly serves 
  to reduce the amount of code needed for HDF related Level 1B
  classes by containing the common values needed.

  With the HDF files, we typically generate the uncertainty through a 
  NoiseModel. Because it is convenient, we include this functionality
  in this class. This is separate from the handling of common values,
  so if needed we can pull this out as a separate mixin. But for now,
  there doesn't seem to be any point to do this.
*******************************************************************/
class Level1bHdf: public Level1b {
public:
  virtual ~Level1bHdf() {}
  virtual int number_spectrometer() const
  { return altitude_.value.rows();}
  virtual DoubleWithUnit latitude(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return latitude_(i);
  }
  virtual DoubleWithUnit longitude(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return longitude_(i);
  }
  virtual DoubleWithUnit solar_zenith(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return solar_zenith_(i);
  }
  virtual DoubleWithUnit solar_azimuth(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return solar_azimuth_(i);
  }
  virtual DoubleWithUnit altitude(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return altitude_(i);
  }
  virtual DoubleWithUnit sounding_zenith(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return sounding_zenith_(i);
  }
  virtual DoubleWithUnit sounding_azimuth(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return sounding_azimuth_(i);
  }
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return stokes_coef_(i, blitz::Range::all());
  }
  virtual DoubleWithUnit relative_velocity(int i) const 
  { return relative_velocity_; }
  virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const
  { ArrayWithUnit<double, 1> res;
    res.value.resize(spectral_coefficient_.value.cols());
    res.value = spectral_coefficient_.value(Spec_index, blitz::Range::all());
    res.units = spectral_coefficient_.units;
    return res;
  }

  virtual Time time(int Spec_index) const { return time_;}
  virtual int64_t sounding_id() const { return hdf_sounding_id_->sounding_id(); }

  // This number is only used to set into the outputted HDF file
  virtual int exposure_index() const 
  { return hdf_sounding_id_->frame_number() + 1; }

  const boost::shared_ptr<HdfFile>& hdf_file_ptr() const {return hfile;}

  boost::shared_ptr<HdfSoundingId> hdf_sounding_id() const 
    {return hdf_sounding_id_;}

  virtual void print(std::ostream& Os) const;

//-----------------------------------------------------------------------
/// Return noise model. Note that this may be a null pointer, which
/// means that we don't have a noise model.
//-----------------------------------------------------------------------

  const boost::shared_ptr<NoiseModel>& noise_model() const 
  { return noise_model_;}

//-----------------------------------------------------------------------
/// Set the noise model
//-----------------------------------------------------------------------
  void noise_model(const boost::shared_ptr<NoiseModel>& Noise_model)
  { noise_model_ = Noise_model;}

  virtual SpectralRange radiance(int Spec_index) const;
protected:
  // Constructor needed when adding new constructor types
  // to inherited classes
  Level1bHdf() { };
  Level1bHdf(const std::string& Fname, 
	     const boost::shared_ptr<HdfSoundingId>& Sounding_id);
  Level1bHdf(const boost::shared_ptr<HdfFile>& Hfile, 
	     const boost::shared_ptr<HdfSoundingId>& Sounding_id);
  std::string file_name;
  boost::shared_ptr<HdfFile> hfile;
  boost::shared_ptr<HdfSoundingId> hdf_sounding_id_;
  ArrayWithUnit<double, 1> altitude_, latitude_, longitude_,
    solar_azimuth_, solar_zenith_, sounding_zenith_, sounding_azimuth_;
  blitz::Array<double, 2> stokes_coef_;
  DoubleWithUnit relative_velocity_;
  ArrayWithUnit<double, 2> spectral_coefficient_;
  Time time_;
  boost::shared_ptr<NoiseModel> noise_model_;

//-----------------------------------------------------------------------
/// Return radiance without the uncertainty filled in. We will fill
/// this in with the noise model, if present.
//-----------------------------------------------------------------------

  virtual SpectralRange radiance_no_uncertainty(int Spec_index) const = 0;
};
}
#endif
