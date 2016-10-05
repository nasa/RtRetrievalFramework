// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "level_1b_average.h"
%}
%base_import(level_1b)
%import "double_with_unit.i"
%import "spectral_range.i"
%import "fp_time.i"
%fp_shared_ptr(FullPhysics::Level1bAverage);

namespace FullPhysics {
class Level1bAverage: public Level1b {
public:
  Level1bAverage(const std::vector<boost::shared_ptr<Level1b> >& Data);
  virtual DoubleWithUnit latitude(int i) const;
  virtual DoubleWithUnit longitude(int i) const;
  virtual DoubleWithUnit sounding_zenith(int i) const;
  virtual DoubleWithUnit sounding_azimuth(int i) const;
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const;
  virtual DoubleWithUnit solar_zenith(int i) const;
  virtual DoubleWithUnit solar_azimuth(int i) const;
  virtual DoubleWithUnit altitude(int i) const;
  virtual DoubleWithUnit relative_velocity(int i) const;
  virtual Time time(int i) const;
  virtual SpectralRange radiance(int Spec_index) const;
};
}
