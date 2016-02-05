// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "level_1b_fts.h"
%}
%base_import(level_1b)
%import "array_with_unit.i"
%import "hdf_file.i"
%import "spectral_bound.i"
%import "double_with_unit.i"
%import "fp_time.i"
%fp_shared_ptr(FullPhysics::Level1bFts);

namespace FullPhysics {
class Level1bFts: public Level1b {
public:
  Level1bFts(const std::string& Runlog,
	     const std::vector<std::string>& Spectra,
	     const ArrayWithUnit<double, 2>& Spectral_range);
  Level1bFts(const std::string& Runlog,
	     const std::vector<std::string>& Spectra,
	     const SpectralBound& Spec_bound);
  Level1bFts(const HdfFile& Hfile,
             const std::vector<std::string>& Band_names,
             const std::string& Radiance_dataset = "/SpectralParameters/modeled_radiance");
  virtual DoubleWithUnit latitude(int i) const;
  virtual DoubleWithUnit longitude(int i) const;
  virtual DoubleWithUnit sounding_zenith(int i) const;
  virtual DoubleWithUnit sounding_azimuth(int i) const;
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const;
  virtual DoubleWithUnit solar_zenith(int i) const;
  virtual DoubleWithUnit solar_azimuth(int i) const;
  virtual DoubleWithUnit altitude(int i) const;
  virtual DoubleWithUnit relative_velocity(int i) const;
  virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const;
  virtual Time time(int Spec_index) const;
  virtual SpectralRange radiance(int Spec_index) const;
};
}
