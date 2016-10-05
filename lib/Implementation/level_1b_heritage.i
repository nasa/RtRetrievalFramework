// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "level_1b_heritage.h"
%}
%base_import(level_1b)
%import "noise_model.i"
%import "double_with_unit.i"
%import "array_with_unit.i"
%import "fp_time.i"
%fp_shared_ptr(FullPhysics::Level1bHeritage);

namespace FullPhysics {
class Level1bHeritage: public Level1b {
public:
  Level1bHeritage(const std::string& Sounding_info_file, 
		  const std::string& Spectrum_file,
		  const boost::shared_ptr<NoiseModel>& Nm);
  virtual DoubleWithUnit latitude(int i) const;
  virtual DoubleWithUnit longitude(int i) const;
  virtual DoubleWithUnit sounding_zenith(int i) const;
  virtual DoubleWithUnit sounding_azimuth(int i) const;
  virtual DoubleWithUnit solar_zenith(int i) const;
  virtual DoubleWithUnit solar_azimuth(int i) const;
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const;
  virtual DoubleWithUnit altitude(int i) const;
  virtual DoubleWithUnit relative_velocity(int i) const;
  virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const;
  virtual Time time(int Spec_index) const;
  virtual SpectralRange radiance(int Spec_index) const;
};
}
