// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "level_1b_oco.h"
%}
%base_import(level_1b_hdf)
%import "hdf_sounding_id.i"
%import "hdf_file.i"
%import "spectral_range.i"
%fp_shared_ptr(FullPhysics::Level1bOco);

namespace FullPhysics {

%feature("notabstract") Level1bOco;

class Level1bOco: public Level1bHdf {
public:
  Level1bOco(const std::string& Fname, 
	     const boost::shared_ptr<HdfSoundingId>& Sounding_id);
  Level1bOco(const boost::shared_ptr<HdfFile>& Hfile, 
	     const boost::shared_ptr<HdfSoundingId>& Sounding_id);
  virtual SpectralRange radiance(int Spec_index) const;
  bool has_spike_eof(int Spec_index) const;
  blitz::Array<double, 1> spike_eof(int Spec_index) const;
  %python_attribute(has_solar_relative_velocity, bool);
  %python_attribute(land_fraction, double);
  %python_attribute(solar_distance, DoubleWithUnit);
  %python_attribute(solar_velocity, DoubleWithUnit);
  %python_attribute(acquisition_mode, std::string);
};
}
