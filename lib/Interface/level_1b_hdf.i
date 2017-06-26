// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "level_1b_hdf.h"
%}

%base_import(level_1b)
%import "noise_model.i"
%import "hdf_sounding_id.i"
%import "hdf_file.i"

%fp_shared_ptr(FullPhysics::Level1bHdf);

namespace FullPhysics {
class Level1bHdf: public Level1b {
public:
  virtual ~Level1bHdf();
  virtual DoubleWithUnit latitude(int i) const;
  virtual DoubleWithUnit longitude(int i) const;
  virtual DoubleWithUnit sounding_zenith(int i) const;
  virtual DoubleWithUnit sounding_azimuth(int i) const;
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const;
  virtual DoubleWithUnit solar_zenith(int i) const;
  virtual DoubleWithUnit solar_azimuth(int i) const;
  virtual DoubleWithUnit altitude(int i) const;
  virtual DoubleWithUnit relative_velocity(int i) const;
  virtual ArrayWithUnit<double, 1> spectral_coefficient(int i) const;
  virtual Time time(int i) const;
  %python_attribute_with_set(noise_model, boost::shared_ptr<NoiseModel>);
protected:
  Level1bHdf();
  Level1bHdf(const std::string& Fname, 
	     const boost::shared_ptr<HdfSoundingId>& Sounding_id);
  Level1bHdf(const boost::shared_ptr<HdfFile>& Hfile, 
	     const boost::shared_ptr<HdfSoundingId>& Sounding_id);
};
}
